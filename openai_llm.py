"""OpenAI call helpers with routing for reasoning vs sampling models."""

from __future__ import annotations

import sys
import threading
from typing import Callable, Optional

_usage_sink: Optional[Callable[[dict], None]] = None
_usage_lock = threading.Lock()


def set_usage_sink(sink: Optional[Callable[[dict], None]]) -> None:
  """Register a callback receiving token usage for every successful LLM call.

  The callback gets a dict with ``model``, ``api``, ``input_tokens``,
  ``cached_input_tokens``, ``output_tokens`` and ``reasoning_tokens``.
  Cached input tokens are a subset of input tokens; reasoning tokens are a
  subset of output tokens. Pass ``None`` to disable tracking.
  """
  global _usage_sink
  with _usage_lock:
    _usage_sink = sink


def _as_int(value) -> int:
  try:
    return int(value or 0)
  except (TypeError, ValueError):
    return 0


def _record_usage(model: str, response, api: str) -> None:
  with _usage_lock:
    sink = _usage_sink
  if sink is None:
    return

  usage = getattr(response, "usage", None)
  if usage is None:
    return

  if api == "responses":
    input_tokens = _as_int(getattr(usage, "input_tokens", None))
    output_tokens = _as_int(getattr(usage, "output_tokens", None))
    input_details = getattr(usage, "input_tokens_details", None)
    output_details = getattr(usage, "output_tokens_details", None)
  else:
    input_tokens = _as_int(getattr(usage, "prompt_tokens", None))
    output_tokens = _as_int(getattr(usage, "completion_tokens", None))
    input_details = getattr(usage, "prompt_tokens_details", None)
    output_details = getattr(usage, "completion_tokens_details", None)

  try:
    sink(
      {
        "model": model,
        "api": api,
        "input_tokens": input_tokens,
        "cached_input_tokens": _as_int(getattr(input_details, "cached_tokens", None)),
        "output_tokens": output_tokens,
        "reasoning_tokens": _as_int(getattr(output_details, "reasoning_tokens", None)),
      }
    )
  except Exception as exc:  # never let accounting break an annotation run
    _print_llm_error("usage sink", model, exc)


def is_reasoning_model(model: str) -> bool:
  m = (model or "").strip().lower()
  if not m:
    return False
  # gpt-5-chat* is a non-reasoning chat model and still accepts sampling params.
  if m.startswith("gpt-5-chat"):
    return False
  return (
    m.startswith("gpt-5")
    or m.startswith("o1")
    or m.startswith("o3")
    or m.startswith("o4")
  )


def resolve_llm_model(raw_model, default: str = "gpt-5.2") -> str:
  if raw_model is None:
    return default
  model = str(raw_model).strip()
  return model or default


def _print_llm_error(stage: str, model: str, exc: BaseException) -> None:
  msg = f"[CellTyper LLM] {stage} failed (model={model}): {type(exc).__name__}: {exc}"
  print(msg, file=sys.stderr, flush=True)


def openai_text_from_messages(
  llm,
  model: str,
  messages,
  *,
  temperature: float = 0.8,
  prompt_cache_key: Optional[str] = None,
  reasoning_effort: str = "medium",
) -> str:
  """Return assistant text, routing around unsupported sampling params.

  Reasoning models (gpt-5.*, o*) reject custom temperature / top_p / penalties.
  For those, call Responses API without temperature. Other models use Chat
  Completions with temperature.

  On API failure, print the original exception message to stderr and try
  fallbacks (still returning "" if all attempts fail).
  """
  model = (model or "").strip()
  if not model:
    return ""

  if is_reasoning_model(model):
    kwargs = {
      "model": model,
      "input": messages,
      "reasoning": {"effort": reasoning_effort},
    }
    if prompt_cache_key:
      kwargs["prompt_cache_key"] = prompt_cache_key
    try:
      response = llm.responses.create(**kwargs)
      _record_usage(model, response, "responses")
      text = getattr(response, "output_text", None) or ""
      if text.strip():
        return text
    except Exception as exc:
      _print_llm_error("responses.create", model, exc)

    # Chat Completions fallback without sampling params (temperature rejected).
    chat_kwargs = {
      "model": model,
      "messages": messages,
    }
    if prompt_cache_key:
      chat_kwargs["extra_body"] = {"prompt_cache_key": prompt_cache_key}
    try:
      completion = llm.chat.completions.create(**chat_kwargs)
      _record_usage(model, completion, "chat")
      text = completion.choices[0].message.content or ""
      if text.strip():
        return text
    except Exception as exc:
      _print_llm_error("chat.completions.create (fallback)", model, exc)

    try:
      response = llm.responses.create(
        model=model,
        input=messages,
        reasoning={"effort": reasoning_effort},
      )
      _record_usage(model, response, "responses")
      return getattr(response, "output_text", None) or ""
    except Exception as exc:
      _print_llm_error("responses.create (retry)", model, exc)
      return ""

  chat_kwargs = {
    "model": model,
    "messages": messages,
    "temperature": temperature,
  }
  if prompt_cache_key:
    chat_kwargs["extra_body"] = {"prompt_cache_key": prompt_cache_key}
  try:
    completion = llm.chat.completions.create(**chat_kwargs)
    _record_usage(model, completion, "chat")
    return completion.choices[0].message.content or ""
  except Exception as exc:
    _print_llm_error("chat.completions.create", model, exc)
    return ""
