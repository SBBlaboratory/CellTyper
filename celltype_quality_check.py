from rapidfuzz import process, fuzz



def quality_check(celltypes, ont, score_cutoff=90):
  terms = []
  for term in ont.terms():
    if term.id.startswith("CL:"):
      terms.append(term.name)
      for syn in term.synonyms:
        terms.append(syn.description)
  
  terms = list(set(terms))
  
  for i in range(len(celltypes)):
    celltype = celltypes[i]
    if type(celltype) != str:
      continue
    celltype = celltype.strip()
    
    flag = True
    
    for term in ont.terms():
      if term.name.lower() == celltype.lower():
        celltypes[i] = term.name
        flag = False
        break
      if any(syn.description.lower() == celltype.lower() for syn in term.synonyms):
        celltypes[i] = term.name
        flag = False
        break
    
    if flag:
      match, score, _ = process.extractOne(celltype.lower(), terms, scorer=fuzz.token_sort_ratio)
      if score >= score_cutoff:
        for term in ont.terms():
          if term.name == match:
            celltypes[i] = term.name
            break
          if any(syn.description == match for syn in term.synonyms):
            celltypes[i] = term.name
            break
      else:
        celltypes[i] = ""
  
  return celltypes



      
            
      
        
        
  
  




