from igraph import Graph
from celltype_quality_check import quality_check
import numpy as np





def construct_marker_graph(state):
  markers = state["markers"]
  species = state["species"]
  ont = state["ont"]
  db1 = state["db1"]
  db2 = state["db2"]
  db3 = state["db3"]
  db8 = state["db8"]
  
  markers = list(markers.split(", "))
  
  graph = Graph(directed=False)
  
  graph.add_vertices(len(markers))
  

  cell_vertex_for_name = {}

  for i in range(len(markers)):
    marker = markers[i]
    graph.vs[i]["name"] = marker
    graph.vs[i]["type"] = "Marker"
    celltypes = []
    celltypes_2 = []
    celltypes_3 = []
    
    if species == "arabidopsis":   
      filtered_db8 = db8.loc[db8["gene"].str.lower() == marker.lower()][:]
      mask = filtered_db8["pct_1"].ge(0.5)
      filtered_db8["reliability"] = ""
      filtered_db8.loc[mask, "reliability"] = np.select(
        [filtered_db8.loc[mask, "pct_diff"].gt(0.5), filtered_db8.loc[mask, "pct_diff"].gt(0.3)],
        ["(high reliable)", "(medium reliable)"],
        default="(low reliable)"
      )
      
      celltypes_2_temp = filtered_db8.loc[filtered_db8["reliability"].ne(""), ["clusterName", "reliability"]]
      
      celltypes_2 = [a + b for a, b in zip(celltypes_2_temp["clusterName"].tolist(), celltypes_2_temp["reliability"].tolist())]
      
      priority = {"high reliable": 3, "medium reliable": 2, "low reliable": 1}

      best = {}
      for s in celltypes_2:
          name, rel = s.rsplit("(", 1)[0].strip(), s.rsplit("(", 1)[1].rstrip(")").strip()
          if name not in best or priority.get(rel, 0) > priority.get(best[name][1], 0):
              best[name] = (s, rel)

      celltypes_2 = [best[name][0] for name in best]
      

      
    else:
      filtered_db1 = db1.loc[(db1["species"].str.lower() == species.lower()) & (db1["marker"].str.lower() == marker.lower())][:]
      celltypes = filtered_db1["cell_name"].unique().tolist()
    
      species_for_db2 = {"human":"hs", "mouse":"mm"}
      filtered_db2 = db2.loc[((species_for_db2[species.lower()] == db2["species"].str.lower()) | (db2["species"].str.lower() == "mm hs")) & (db2["official gene symbol"].str.lower() == marker.lower())][:]
      celltypes_2 = filtered_db2["cell type"].unique().tolist()
    
    
    if species != "arabidopsis":
      species_for_db3 = {"human":"Homo sapiens", "mouse":"mouse"}
      db3["extra"] = db3["gene_symbol"].str.lower().str.split(r",\s*")
      db3["extra"] = db3["extra"].apply(lambda l: l if isinstance(l, list) else [])
    
        
      filtered_db3 = db3.loc[(db3["species"].str.lower() == species_for_db3[species.lower()].lower()) & (db3["extra"].apply(lambda l: marker.lower() in l))][:]
      celltypes_3 = filtered_db3["cell_type"].unique().tolist()

    celltypes.extend(celltypes_2)
    celltypes.extend(celltypes_3)
    celltypes = list(set(celltypes))
    
    
    if species == "human":
      celltypes = list(set(quality_check(celltypes, ont)))
    
    for j in range(len(celltypes)):
      celltype = celltypes[j]
      if type(celltype) != str or celltype == "":
        continue
      if celltype not in cell_vertex_for_name:
        graph.add_vertices(1)
        vid = graph.vcount() - 1
        graph.vs[vid]["name"] = celltype
        graph.vs[vid]["type"] = "Cell"
        cell_vertex_for_name[celltype] = vid
      graph.add_edge(i, cell_vertex_for_name[celltype])
  
  isolated = [v.index for v in graph.vs if v.degree() == 0]
  if isolated:
    graph.delete_vertices(isolated)
  pbar = state.get("pbar")
  if pbar is not None:
    try:
      pbar.update(20)
    except Exception:
      pass
  return {**state, "marker_cell_graph":graph, "pbar":pbar}





