#' @title umi_sim_graph
#'
#' @description build UMI graph based on their sequence similarity
#' @details The function calculates the similarity of two UMIs by needleman score. Their isoform similarity will also be considered if given.
#' Based on the similarity a sNN graph will be built.
#' @param umi The umi sequence for each read
#' @param iso The isoform sequence for each read
#' @param sim_thresh The minimum threshold of the similarity. Only UMIs with higher similarity will be connected with an edge.
#' @param iso_thresh The maximum threshold that two isoforms can be different in overlapping bases.
#' @param split The character to split each exon in the isoform representation
#' @param sep The character to seperate the exon start and end in the isoform representation
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph simplify
#' @importFrom igraph set_vertex_attr
#'
umi_sim_graph <- function(umi,iso = NULL,sim_thresh = 7,
                              iso_thresh = 200,split = "|",sep = ","){
  if(is.null(iso)){
    iso = rep("N",length(umi))
  }
  else if(length(umi) != length(iso)){
    warning("The length of umi and isoforms don't correspond,
                the information of isoforms won't be used!")
    iso = rep("N",length(umi))
  }

  iso_umi_table = as.data.frame(cbind(iso,umi))
  colnames(iso_umi_table) = c("iso","umi")
  iso_umi_table$id = 1:nrow(iso_umi_table)
  iso_umi_table = iso_umi_table %>%
                  group_by(iso,umi) %>%
                  summarise(count = n(),id = list(id),.groups = "drop")
  #return(iso_umi_table)
  umi_ns = umi_graph_table(iso_umi_table$umi,iso_umi_table$iso,iso_umi_table$count,
                           sim_thresh,iso_thresh,split,sep);

  umi_ns = as.data.frame(do.call(rbind,umi_ns))
  colnames(umi_ns) <- c("node1","node2","ns","count")

  umi_ns = umi_ns[umi_ns$ns > 0,]
  umi_ns$node1 = umi_ns$node1+1
  umi_ns$node2 = umi_ns$node2+1
  umi_ns = umi_ns %>% filter(count > 0) %>% mutate(weight = ns)
  # return(umi_ns)

  graph = graph_from_data_frame(umi_ns,directed = FALSE,
                                vertices = 1:nrow(iso_umi_table))
  graph = igraph::simplify(graph,remove.loops = TRUE,edge.attr.comb = "first")
  graph = graph %>%
    set_vertex_attr("count", value = iso_umi_table$count)

  umi_corres_size = sapply(iso_umi_table$id,length)
  umi_corres = rep(1:nrow(iso_umi_table),umi_corres_size)
  umi_corres = as.data.frame(cbind(do.call(c,iso_umi_table$id),umi_corres))
  colnames(umi_corres) = c("id","pair_id")
  return(list(graph,umi_corres))
}


#' @title createSNN
#'
#' @description count the number of shared neighbors between two connected nodes in a UMI graph
#' @details The function receives the edge table as input, count the shared neighbors between two connected nodes
#' @param dis The edge table from a UMI graph with two columns, each col represents a node.
#' @param count The weight for each node, if NULL, all weight will be set as 1.
#' @param self For the node with self loop in the graph, decide if themselves should be considered in the neibor counting.
#' @param iso_thresh The maximum threshold that two isoforms can be different in overlapping bases.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
createSNN <- function(dis,count = NULL,self = FALSE){
  dis <- as.matrix(dis)
  dis <- rbind(dis,dis[,c(2,1)])

  dis <- as.data.frame(dis)
  colnames(dis) <- c("node1","node2")
  dis$node1 = as.character(dis$node1)
  dis$node2 = as.character(dis$node2)

  if(!self){
    dis = dis %>% filter(node1 != node2)
  }
  if(is.null(count)){
    count = rep(1,length(unique(dis$node1)))
    names(count) = unique(dis$node1)
  }

  neighbors <- as.data.frame(dis %>% group_by(node1) %>%
                               summarise(neighbor = list(node2)))

  index = neighbors$node1
  neighbor = neighbors$neighbor
  count = count[index]

  snn <- shareNeighbor(index,neighbor,count)

  return(snn)
}

#' @title SNN_graph
#'
#' @description add neighbor sharing as an edge attribute in a graph
#' @details The function calculates the neighbor sharing for a graph, and add the count of shared neighbors as an edge attribute in the graph
#' @param graph an igraph object
#' @importFrom igraph as_long_data_frame
#' @importFrom igraph edge_attr
#' @importFrom igraph vertex_attr
SNN_graph <- function(graph){
    vertex = vertex_attr(graph)$name
    count = vertex_attr(graph)$count
    names(count) = vertex

    edge_table = as_long_data_frame(graph)[,c("from","to")]
    edge_table = cbind(vertex[edge_table$from],vertex[edge_table$to])

    snn = createSNN(dis = edge_table,count = count,self = FALSE)
    snn <- snn[snn$node1 <= snn$node2,]

    index <- paste(snn$node1,snn$node2,sep = "|")
    edge_attr(graph, "share", index = index) <- as.numeric(snn$share)
    edge_attr(graph)$share[is.na(edge_attr(graph)$share)] = 1

    return(graph)
}

#' @title louvain_iter_stack
#'
#' @description do iterative louvain clustering for nodes in a graph
#' @details The function does iterative louvain clustering for an igraph object. Clusters will be
#' split as sub-graphs and stores in a list.
#' @param graph an igraph object.
#' @param weight the name of the edge attribute to use as weight in louvain clustering.
#' @param resolution the resolution for the louvain clustering.
#' @param alpha the fold change between the number of edges and number of nodes. Lower alpha requires
#' highlier connected graph.
#' @param sim_thresh the minimum threshold for the needleman similarity for a graph to stop clustering
#' @importFrom igraph decompose
#' @importFrom igraph graph.mincut
#' @importFrom igraph edge_attr
#' @importFrom igraph cluster_leiden
#' @importFrom igraph delete_edges
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph crossing
louvain_iter_stack <- function(graph,weight = "weight",
                         resolution = 1,alpha = 2, sim_thresh = 6){
  in_graph_list = decompose(graph)
  out_graph_list = list()

  while(length(in_graph_list) > 0){
    temp_graph = in_graph_list[[1]]
    in_graph_list[[1]] = NULL

    if(length(V(temp_graph)) == 1){
      out_graph_list = append(out_graph_list,list(temp_graph))
    }

  else{
    min_cut <- graph.mincut(temp_graph)
    if(min_cut >= length(V(temp_graph))/alpha |
       min(edge_attr(temp_graph)[[weight]]) > sim_thresh){
      out_graph_list = append(out_graph_list,list(temp_graph))
    }
    else{
      #temp_graph = SNN_graph(temp_graph)

      #edge_weight <- edge_attr(temp_graph)[[weight]]
      #edge_share <- edge_attr(temp_graph)[["share"]]
      #edge_share[edge_weight > sim_thresh & edge_share == 0] <- 1

      #edge_attr(temp_graph,"weight",index = E(temp_graph)) = edge_weight*sqrt(edge_share)
      cluster <- cluster_leiden(temp_graph,resolution = resolution)

      graph_cut <- delete_edges(temp_graph, E(temp_graph)[crossing(cluster,temp_graph)])

      sub_graphs = decompose(graph_cut)

      if(length(sub_graphs) == 1){
        out_graph_list = append(out_graph_list,sub_graphs)
      }
      else{
        in_graph_list = append(in_graph_list, sub_graphs)
      }
    }
  }
  }
  return(out_graph_list)
}
#' @title graph_to_cluster
#'
#' @description extract the cluster information for nodes from a list of sub-graphs.
#' @details each sub-graph in the list represents a cluster, this function assigns the cluster id
#' for nodes within those sub-graphs.
#' @param graph_list a list of igraph object, each object represents a cluster.
#' @importFrom igraph vertex_attr
graph_to_cluster <- function(graph_list){
    out <- lapply(1:length(graph_list),function(i){
        cluster = graph_list[[i]]
        cluster_table <- cbind(vertex_attr(cluster)$name,i)
    })
    out <- as.data.frame(do.call(rbind,out))
    colnames(out) <- c("id","cluster")
    out$id = as.numeric(out$id)
    out$cluster = as.numeric(out$cluster)
    out <- out[order(out[,"id"]),]

    return(out[,"cluster"])
}

#' @title umi_cluster
#'
#' @description apply iterative louvain clustering for UMIs.
#' @details apply iterative louvain clustering for UMIs based on the needleman score between each UMI. Isoform
#' information is optional to be used to help distinguish different UMIs.
#' @param umi a vector of UMIs.
#' @param iso a vector of isoform sequences.
#' @param sim_thresh the minimum threshold of needleman score to connect to different UMIs.
#' @importFrom dplyr left_join
umi_cluster <- function(umi,iso = NULL,thresh = NULL){
  if(is.null(thresh)){
    thresh = round(length(umi)/2+1)
  }
  graph_corres <- umi_sim_graph(umi,iso = iso,sim_thresh = thresh)
  graph = graph_corres[[1]]
  umi_corres = graph_corres[[2]]

  graph_cluster = louvain_iter_stack(graph = graph,alpha = 2,sim_thresh = thresh+2)
  cluster = graph_to_cluster(graph_cluster)
  cluster = as.data.frame(cbind(1:length(cluster),cluster))
  colnames(cluster) = c("pair_id","cluster")

  cluster_expand = left_join(umi_corres,cluster,by = "pair_id")
  cluster_expand = cluster_expand[order(cluster_expand$id),]
  return(cluster_expand$cluster)
}
