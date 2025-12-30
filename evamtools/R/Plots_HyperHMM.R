#FUNCIONES PARA PLOTEAR HYPERHMM

#Código para generar plots extraído de:
# plot_HyperHMM_bubbles ->  hyperhmm::plot_bubbles
# plot_HyperHMM_pfg ->  hyperhmm::plot_pfg
# plot_HyperHMM_hypercube ->  hyperhmm::plot_hypercube
# plot_HyperHMM_standard ->  hyperhmm::plot_standard
# plot_HyperHMM_hypercube_flux ->  hyperhmm::plot_hypercube_flux

#BUBBLE PLOT
plot_HyperHMM_bubbles <- function(fitted.obj, # output data structure. either just a matrix of probabilities (formatted == F) or a dataframe output from the HyperHMM wrapper with means and sds (formatted == T)
                                    labels = NULL, # labels for feature names
                                    formatted = FALSE) # dataframe formatted or not? see above
{  
    message("Building bubble plot")

  bp = fitted.obj$stats
  # if we've just got a matrix of probabilities, pull it into long form
  if(formatted == FALSE) {
    bp.df = data.frame()
    for(i in 1:nrow(bp)) {
      for(j in 1:ncol(bp)) {
        bp.df = rbind(bp.df, data.frame(order=j, feature=i, prob=bp[i,j]))
      }
    }
  } else { bp.df = bp; bp.df$prob = bp.df$mean  }
  ## plot bubbles
  if(is.null(labels)) {
    g.1 = ggplot2::ggplot(bp.df, ggplot2::aes(x=order, y=feature)) +
      ggplot2::geom_point(ggplot2::aes(size=prob), colour="#31b786ff") +
      ggplot2::scale_x_continuous(breaks=1:max(bp.df$order)) +
      ggplot2::scale_y_continuous(breaks=1:max(bp.df$feature)) +
      ggplot2::theme_classic() + ggplot2::theme(legend.position = "none")
  } else {
    g.1 = ggplot2::ggplot(bp.df, ggplot2::aes(x=order, y=feature)) +
      ggplot2::geom_point(ggplot2::aes(size=prob), colour="#000000ff") +
      ggplot2::scale_y_continuous(breaks=length(labels):1, labels=labels) +
      ggplot2::theme_classic() + ggplot2::theme(legend.position = "none")
  }
  return(g.1)
}

#-----------PFG PLOT
plot_HyperHMM_pfg <- function(fitted.obj,  # list of transitions between states
                    pfg.layout = "matrix",  # graph layout
                    curvature = 1)   # geometric parameter for edge curviness
{

  translist = fitted.obj$viz
  message("Building PFG")
  ## PFG
  # get the set of pairwise first-next feature labels throughout each pathway
  edges=data.frame()
  message("Reading")
  zeroes = strsplit(translist[1], split=" ")[[1]][1]
  for(i in 1:1000) {
    #print(i)
    # get source and destination states
    src = strsplit(translist[i], " ")[[1]][1]
    dest = strsplit(translist[i], " ")[[1]][2]
    srcn = as.numeric(strsplit(src, "")[[1]])
    destn = as.numeric(strsplit(dest, "")[[1]])
    # identify changed feature
    change = which(srcn-destn != 0)
    # add this and last change to list of ordered pairs
    if(length(change) == 1) {
      if(src != zeroes) {
        edges = rbind(edges, data.frame(src=lastchange, dest=change))
      } else {
        edges = rbind(edges, data.frame(src=0, dest=change))
      }
      lastchange = change
    }
  }
  # get unique edges in this adjacency list and counts
  uedges = unique(edges)
  ucount = 0*uedges$src
  for(i in 1:nrow(uedges)) {
    ucount[i] = length(which(edges$src == uedges$src[i] & edges$dest == uedges$dest[i]))
  }

  message("Building adj mat")
  # construct graph from these edges
  uedges = rbind(uedges, uedges[1,])
  ucount[length(ucount)+1] = 1
  g = igraph::graph.data.frame(uedges)
  igraph::E(g)$weight = ucount
  maxw = max(ucount)
  sumw = sum(ucount)

  # sort nodes to a canonical order -- helps when we're comparing different cases, especially with the matrix layout
  s <- gtools::mixedsort(names(igraph::V(g)))
  new.g = igraph::permute(g, match(igraph::V(g)$name, s))
  g = new.g
  igraph::V(g)$name[1] = "-"

  message("Building plot")
  # plot PFG
  if(pfg.layout == "tree") {
    g.3 = ggraph::ggraph(g, layout="tree") +
      ggraph::geom_edge_bend(ggplot2::aes(edge_width=exp(weight/sumw), edge_alpha = weight/sumw),
                             strength=curvature,  arrow=ggplot2::arrow()) +
      ggraph::geom_node_point() +
      ggraph::geom_node_label(ggplot2::aes(label=name), nudge_x = 0.05, nudge_y=-0.05) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none")
  } else if(pfg.layout == "matrix") {
    g.3 = ggraph::ggraph(g, layout="matrix") +
      ggraph::geom_edge_bend(ggplot2::aes(edge_width=exp(weight/sumw), edge_alpha = weight/sumw, color="#000000cc"),
                             strength=curvature,  arrow=ggplot2::arrow()) +
      ggraph::geom_node_point(ggplot2::aes(color=name)) +
# ggraph::geom_node_label(ggplot2::aes(label=name, fill=name))  # fondo de etiqueta según nombre
      ggraph::geom_node_label(ggplot2::aes(label=name, fill=name), nudge_x = 0.05, nudge_y=-0.05) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none")
  } else {
    g.3 = ggraph::ggraph(g) +
      ggraph::geom_edge_bend(ggplot2::aes(edge_width=exp(weight/sumw), edge_alpha = weight/sumw),
                             strength=curvature,  arrow=ggplot2::arrow()) +
      ggraph::geom_node_point() +
      ggraph::geom_node_label(ggplot2::aes(label=name), nudge_x = 0.05, nudge_y=-0.05) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none")
  }
  return(g.3)
}

#HYPERCUBE PLOT
plot_HyperHMM_hypercube <- function(fitted.obj,               # including set of transitions
                           use.width = TRUE,           # use line width to display edge weights?
                           duplicate.offset = 0.,   # vertical offset for nodes in identical positions
                           lab.size = 3,            # size for edge labels
                           p.size = 1,              # point size
                           node.labels = TRUE,         # node labels, yes or no?
                           seg.labels = TRUE,          # line segment labels?
                           threshold = 0,           # ignore edges under a threshold in the hypercube plot
                           break.redundancy = 0,    # itself redundant now?
                           rotate.phi = FALSE)          # rotate states out of the page (in case of trajectories bunched up near the top/bottom)
{
  translist = fitted.obj$viz
  message("Building hypercube plot")

  ## hypercube
  # get unique set of transitions and associated counts
  l = unique(translist)
  counts = rep(0, length(l))
  for(i in 1:length(l)) {
    set = which(translist == l[i])
    counts[i] = length(set)
  }
  #  l = l[counts > threshold]

  # split into lists of source and destination nodes
  srcs = dests = list()
  n = 1
  for(line in l) {
    s = strsplit(line, " ")
    srcs[[n]] = s[[1]][1]
    dests[[n]] = s[[1]][2]
    n = n + 1
  }
  # set string length and 0^L string
  len = nchar(srcs[[1]])
  zero = paste(rep("0", len), collapse="")

  # produce useful vectors
  srcs = unlist(srcs)
  dests = unlist(dests)

  all.nodes = apply(expand.grid(rep(list(0:1),len)), 1, paste, collapse="")
  nodes = all.nodes
  nnodes = length(nodes)

  # produce list storing where incoming edges to each node come from
  ins = list()
  for(node in nodes) {
    refs = which(dests == node)
    refcodes = srcs[refs]
    ins[[length(ins)+1]] = which(nodes %in% refcodes)
  }

  # produce hypercube visualisation

  message("Calculating embedding")
  # spherical polars: r, theta, phi
  # r = 1 everywhere
  rs = rep(1, nnodes)
  # theta is just set by number of 1s in a string

  thetas = unlist(lapply(nodes, function(s) { return(stringr::str_count(s, "0")*3.14159/len) }) )

  # initialise phis
  phis = rep(-1, nnodes)
  zero.counts = unlist(lapply(nodes, function(s) { return(stringr::str_count(s, "0")) }) )
  for(zeroes in 0:len) {
    refs = which(zero.counts == zeroes)
    these.phis = (0:(length(refs)-1))/length(refs)*3.14159
    phis[refs] = these.phis
  }

  # dataframes for spherical and cartesian coordinates
  spcoords = data.frame(r = rs, theta = thetas, phi = phis, label = nodes)
  # rotate phi values if required
  if(rotate.phi == TRUE) {
    spcoords$phi = spcoords$phi + 3.14159/2
  }

  coords = data.frame(x = spcoords$r*cos(spcoords$phi)*sin(spcoords$theta),
                      y = spcoords$r*sin(spcoords$phi)*sin(spcoords$theta),
                      z = spcoords$r*cos(spcoords$theta), label = spcoords$label)

  # dataframe for line segments in cartesian coords
  segments = data.frame()
  seglabels = data.frame()
  safe.nodes = rep(F, nrow(coords))
  for(i in 1:length(srcs)) {
    if(counts[i] > threshold) {
      src = which(coords$label == srcs[i])
      dest = which(coords$label == dests[i])
      safe.nodes[src] = safe.nodes[dest] = T
      label = paste(c("+", which(unlist(strsplit(srcs[i], split="")) !=unlist(strsplit(dests[i], split="")))), collapse="")
      segment = data.frame(src.x = coords[src,]$x, src.y = coords[src,]$y, src.z = coords[src,]$z, dest.x = coords[dest,]$x, dest.y = coords[dest,]$y, dest.z = coords[dest,]$z, count=counts[i])
      segments = rbind(segments, segment)
      seglabels = rbind(seglabels, data.frame(x = (segment$src.x + segment$dest.x)/2, y = (segment$src.y + segment$dest.y)/2,  z = (segment$src.z + segment$dest.z)/2, label=label))
    }
  }

  base = data.frame(src.x=0,src.z=0,dest.x=0,dest.z=0,count=0)

  message("Make plot")
  # plot
  if(use.width) {
    cube.plot = ggplot2::ggplot() +
      ggplot2::geom_segment(data=segments, ggplot2::aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x, size=count/max(count)), color="#6f0057ff") +
      ggplot2::scale_size_identity() +
      ggplot2::geom_point(data = coords[safe.nodes == T,], ggplot2::aes(x=z, y=x), size=p.size, color="#b107bdff") +
      ggplot2::xlim(-1,1) + ggplot2::ylim(-1,1) +
      ggplot2::theme_void() + ggplot2::theme(legend.position="none")
    if(node.labels == TRUE) {
      cube.plot = cube.plot +
        ggrepel::geom_text_repel(data = coords[safe.nodes == T,], ggplot2::aes(x=z,y=x,label=label))
    }
    if(seg.labels == TRUE) {
      cube.plot = cube.plot +
        ggplot2::geom_text(data=seglabels, ggplot2::aes(x=z,y=x,label=label), color="#136d4cff", size=lab.size)
    }
  } else {
    cube.plot = ggplot2::ggplot() +
      ggplot2::geom_segment(data=segments, ggplot2::aes(x=src.z, y=src.x, xend=dest.z, yend=dest.x), color="#CCCCCC") +
      ggplot2::geom_point(data = coords[safe.nodes == T,], ggplot2::aes(x=z, y=x), size=p.size, color="#ab0505ff") +
      ggplot2::xlim(-1,1) + ggplot2::ylim(-1,1) +
      ggplot2::theme_void() + ggplot2::theme(legend.position="none")
    if(node.labels == TRUE) {
      cube.plot = cube.plot +
        ggrepel::geom_text_repel(data = coords[safe.nodes == T,], ggplot2::aes(x=z,y=x,label=label))
    }
    if(seg.labels == TRUE) {
      cube.plot = cube.plot +
        ggplot2::geom_text(data=seglabels, ggplot2::aes(x=z,y=x,label=label), color="#888888", size=lab.size)
    }
  }
  return(cube.plot)
}

#STANDARD PLOT
plot_HyperHMM_standard <- function(fitted.obj, 
                                    legacy=FALSE, 
                                    label="")
{
    plot_bubs = plot_HyperHMM_bubbles(fitted.obj, formatted=TRUE) + ggplot2::ggtitle(label)
    plot_flux = plot_HyperHMM_hypercube_flux(fitted.obj, thresh = 0.02) +
        ggplot2::theme(legend.position = "none")
    plot_diag = plot_HyperHMM_pfg(fitted.obj, pfg.layout="matrix")
    plot_standard = ggpubr::ggarrange(plot_flux, plot_bubs, plot_diag, nrow=1)

  return(plot_standard)
}

#-----HYPERCUBE FLUX PLOT
#Función necesaria para el plot: hypercube flux
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

plot_HyperHMM_hypercube_flux <- function(fitted.obj, 
                                        thresh = 0.05, 
                                        node.labels = TRUE, 
                                        use.probability = FALSE)
{
    my.post = fitted.obj
  ### produce hypercube subgraph
  bigL = my.post$L
  if(use.probability == TRUE) {
    trans.p = my.post$transitions[my.post$transitions$Probability > thresh & my.post$transitions$Bootstrap == 0,]
  } else {
    trans.p = my.post$transitions[my.post$transitions$Flux > thresh & my.post$transitions$Bootstrap == 0,]
  }
  trans.g = igraph::graph_from_data_frame(trans.p[,2:ncol(trans.p)])
  bs = unlist(lapply(as.numeric(igraph::V(trans.g)$name), DecToBin, len=bigL))
  igraph::V(trans.g)$binname = igraph::V(trans.g)$name
  layers = stringr::str_count(bs, "1")
  if(use.probability == TRUE) {
    this.plot = ggraph::ggraph(trans.g, layout="sugiyama", layers=layers) +
      ggraph::geom_edge_link(ggplot2::aes(edge_width=Probability, edge_alpha=Probability)) +
      ggraph::scale_edge_width(limits=c(0,NA)) +
      ggraph::scale_edge_alpha(limits=c(0,0.5)) +
      ggplot2::theme_void()
  } else {
    this.plot = ggraph::ggraph(trans.g, layout="sugiyama", layers=layers) +
      ggraph::geom_edge_link(ggplot2::aes(edge_width=Flux, edge_alpha=Flux), color="#079760cc") +
      ggraph::scale_edge_width(limits=c(0,NA)) +
      ggraph::scale_edge_alpha(limits=c(0,0.5)) +
      ggplot2::theme_void()
  }
  if(node.labels == TRUE) {
    this.plot = this.plot + ggraph::geom_node_text(ggplot2::aes(label=binname),angle=45,size=2)
  }
  return(this.plot)
}
