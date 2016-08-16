using NetworkViz
using LightGraphs
using Distributions
using Escher

function netgraph(g)
  nvg = nv(g)
  dists = map(x->Normal(x, 0.25), 0.25:0.25:0.75)
  norm = pdf(dists[1], 0.25)
  c = Colors.Color[Colors.RGBA(map(x->pdf(x, i/nvg)/norm, dists)...) for i in 1:nvg]
  e = EdgeProperty("#000000",0.7)
  drawGraph(g, node=NodeProperty(c, 0.1, 1), edge=e, z=1)
end

function main(window)
  push!(window.assets, "widgets")
  push!(window.assets, "tex")
  push!(window.assets,("ThreeJS","threejs"))
  push!(window.assets, "layout2")
  g = DiGraph(100, 300)

  # create the tabs
  tabbar = tabs([
       hbox("Description"),
       hbox("Interactive"),
       hbox("Necklace graph"),
  ])

  # create the pages
  tabcontent = pages([
    title(3, "web component all the things"),
    netgraph(g),
    Escher.image("http://imgur.com/qnPxe3a.png")
  ])

  # connect the tabs to pages
  # returns a pair of "connected" tab set and pages
  t, p = wire(tabbar, tabcontent, :tab_channel, :selected)

  # stack them on top of each other
  vbox(t, p)
end
