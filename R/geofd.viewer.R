".geofd.viewer"<-function(okfd.res){

  # Loading required libraries
  require(tkrplot)

  tt <- tktoplevel()

  coordsFrame <- tkframe(tt,relief="groove",borderwidth=2,width=486)
  predictionsFrame <- tkframe(tt,relief="flat",borderwidth=2,width=486)
  predictionsFrame1 <- tkframe(tt,relief="flat",borderwidth=2,width=486)

  tkwm.title(tt,"geofd results viewer")
  parPlotSize <- c()
  usrCoords <- c()
  index <- 0
  ylim <- c(min(okfd.res[["krig.new.data"]]),max(okfd.res[["krig.new.data"]]))

  # Plot a coordinate map for the specified set of points
  plotCoords <- function()
  {
    params <- par(bg="white")
    plot(okfd.res[["new.coords"]][,1], okfd.res[["new.coords"]][,2], xlab="Xcoord", ylab="Ycoord")
    if(index>=0){
      points(okfd.res[["new.coords"]][index,1], okfd.res[["new.coords"]][index,2],col="red",pch=20,cex=2)
    }
    parPlotSize <<- par("plt")
    usrCoords <<- par("usr")
    par(params)
  }

  # Plot a graph for the specified predicted curve
  plotPredictions <- function()
  {
    params <- par(bg="white")
    if(index>0){
    plot(okfd.res[["argvals"]], okfd.res[["krig.new.data"]][,index], col=1, lwd=1, type="l", lty=1, main="Predictions", xlab="argnames", ylab="variable", ylim=ylim )
    }
  }

  # Handle the click event on the coords plot
  # Determines in plot coordinates the closest point to the clicked one
  # and replot each graph, which variates with the index argument
  OnLeftClick <- function(x,y)
  {
    xClick <- x
    yClick <- y

    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",coordsPlot)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",coordsPlot)))

    xMin <- parPlotSize[1] * width
    xMax <- parPlotSize[2] * width
    yMin <- parPlotSize[3] * height
    yMax <- parPlotSize[4] * height

    rangeX <- usrCoords[2] - usrCoords[1]
    rangeY <- usrCoords[4] - usrCoords[3]

    xClick <- as.numeric(xClick)+0.5
    yClick <- as.numeric(yClick)+0.5
    yClick <- height - yClick

    xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
    yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)

    closest <- closestPoint(xPlotCoord, yPlotCoord)

    # If it is the first time this function is run
    if(index==0){
      # Destroy the message frame and activates the predictions plot frame
      tkdestroy(predictionsFrame)
      tkpack(predictionsFrame1,side='left')
    }

    # Retrieve and save the index of the closest point to the clicked one
    index <<- closest$index

    # Replot coords and predicions plots
    tkrreplot(coordsPlot)
    tkrreplot(predictionsPlot)
  }

  # Returns a data frame with coordinate index, X coord and Y coord
  closestPoint <- function(x, y){
  #  result <- c()
    index <- 0
    mindist <- 9999999999999999
    for(coordnum in seq(1, dim(okfd.res$new.coords)[1])){
      pointdist <- dist(rbind(c(x, y), okfd.res$new.coords[coordnum,]))
      if(pointdist<mindist){
        index <- coordnum
        mindist <- pointdist
      }
  #    result <- rbind(result, c(coordnum, okfd.res$new.coords[coordnum,], pointdist) )
    }
  #  print(result)
    data.frame(index=index, x=as.numeric(okfd.res$new.coords[index,][1]), y=as.numeric(okfd.res$new.coords[index,][2]))
  }

  # Prepare coordinates plot frame
  coordsPlot <- tkrplot(coordsFrame,fun=plotCoords,hscale=1,vscale=1)
  tkgrid(coordsPlot)

  # Prepare message frame
  tkgrid(tklabel(predictionsFrame,text="Please click near a point of your interest on the left frame"))

  # Prepare predictions plot frame
  predictionsPlot <- tkrplot(predictionsFrame1, fun=plotPredictions, hscale=1, vscale=1)
  tkgrid(predictionsPlot)

  # Binding the click event 
  tkbind(coordsPlot, "<Button-1>", OnLeftClick)
  tkconfigure(coordsPlot,cursor="hand2")

  # Packaging the frame
  tkpack(coordsFrame, side='left')
  tkpack(predictionsFrame, side='left')

}
