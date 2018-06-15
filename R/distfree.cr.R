#get the cross points of lines a=slop, b=intercept
linescrosspoint<-function(a, b, a2=NA, b2=NA){
  if (missing(b2)){
    b1=b[-NROW(b)]
    b2=b[-1]
  } else b1=b;
  if (missing (a2)){
    a1=a[-NROW(a)]
    a2=a[-1]
  } else {a1=a}
  
  x1=-(b1-b2)/(a1-a2)
  y1=-(a2*b1-a1*b2)/(a1-a2)
  na.exclude(data.frame(x=x1, y=y1))
}

AngleToSlope<-function(x) {
  if (cos(x)==0) res=1e50 else res=sin(x)/cos(x)
  res
}

linescrosspoint.plus<-function(lines){
  p0=(colMeans(linescrosspoint(a=c(lines[,1], lines[1:2,1]), b=c(lines[,2], lines[1:2,2]))))
  p0=data.frame(x=p0[1], y=p0[2])
  cross.distance=(lines$b*p0$x-p0$y+lines$incpt)/sqrt(lines$b^2+1)
  
  idx=which(abs(cross.distance)==min(abs(cross.distance)))
  line.idx=idx
  
  p1=linescrosspoint(c(lines$b[idx], AngleToSlope(pi/2+lines$angle[idx])),
                     c(lines$incpt[idx], -AngleToSlope(pi/2+lines$angle[idx])*p0$x+p0$y))
  while (1==1) {
    sub.line.idx=(1:((NROW(lines))/4))+line.idx[NROW(line.idx)]
    sub.line.idx[which(sub.line.idx>NROW(lines))]=sub.line.idx[which(sub.line.idx>NROW(lines))]-NROW(lines)
    sub.line=lines[sub.line.idx,]
    cur.line=lines[line.idx[NROW(line.idx)],]
    xy=linescrosspoint(sub.line$b, sub.line$incpt, cur.line$b, cur.line$incpt)
    cross.distance=(xy$x-p1$x)^2+(xy$y-p1$y)^2
    idx=which(abs(cross.distance)==min(abs(cross.distance)))
    idx=idx[NROW(idx)]
    p1=xy[idx,]
    line.idx=c(line.idx, sub.line.idx[idx])
    if (line.idx[NROW(line.idx)]==line.idx[1]) break
  }
  
  linescrosspoint(lines$b[line.idx], lines$incpt[line.idx])
}



#get range quantile for scatter of x, y
xy.quantile<-function(x, y, alpha=0.05, half.nknots=90) { 
  
  probs=range(c(alpha/2, 1-alpha/2))
  anglelist=seq(-pi/2, pi/2-pi*(0.5/half.nknots), length.out=half.nknots+1)
  res=c(1)[0]
  off.points=rep(0, NROW(x))
  for (angle in anglelist) {
    if (angle%%(pi/2)!=0) {
      a=tan(angle)
      b=-1
      d=-1*(a*x+b*y)/sqrt(a^2+b^2)
      crt=quantile(x=d, probs=probs)
      idx=which((d<crt[1])|(d>crt[2]), arr.ind=T)
      off.points[idx]=1
      res=rbind(res, c(a, crt*sqrt(1+a^2), angle, 1))
    }
  }
  
  
  lines=data.frame(b=c(res[,1],res[,1]), incpt=c(res[,2],res[,3]), angle=c(res[,4], res[,4]+pi)%%(2*pi), c=c(res[,5], res[,5]))
#   p=linescrosspoint(a=c(lines[,1], lines[1:2,1]), b=c(lines[,2], lines[1:2,2]))
  p=linescrosspoint.plus(lines)
  
  
  list(polygon=p, data=data.frame(x=x, y=y, pip=1-off.points), alpha.realized=sum(off.points)/NROW(x))
}

plot.distfree.cr<-function(x, show.points=T, ...) {
  obj=x;
  xlim=range(c(obj$data[,1], obj$polygon$x))
  ylim=range(c(obj$data[,2], obj$polygon$y))
  plot(NA, type="n", xlim=xlim, ylim=ylim, xlab=obj$xlab, ylab=obj$ylab)
  if (show.points) {
    points(obj$data, col=obj$col[2:3][2-obj$data$pip])
  }
  polygon(obj$polygon, border=obj$col[1], lwd=2.0)
  legend("topright",lty=c("solid", "blank", "blank"), bty="n",
         legend=c(paste("Confidence region (p=", 1-x$alpha.realized, ")", sep=""),
                  "PIP (points in polygon)", "non-PIP"),
         pch=c(NA, 1, 1),
         col=c("red", "black", "gray")
  )
}

distfree.cr<-function(x, y, alpha=0.05, alpha.min.diff=0.5/NROW(x), nknots=40,
                      xlab = deparse(substitute(x)), 
                      ylab = deparse(substitute(y)),
                      col=c("red", "black", "gray"), draw=T) {

  if (missing(y)) {
    if (NCOL(x) == 2) {
      if (missing(xlab)) xlab <- colnames(x)[1]
      if (missing(ylab)) ylab <- colnames(x)[2]
      y <- x[, 2]
      x <- x[, 1]
    } else stop("x and y must be vectors, or x must be a 2 column matrix")
  } else if (!(is.vector(x) && is.vector(y) && length(x) == length(y))) stop("x and y must be vectors of the same length")
  col=rep(col, 3)[1:3]
  
#   if ((NROW(x)*alpha)<2) {
#     xy.quant=(list(polygon=border.points.idx(data.frame(x=x, y=y)), data=data.frame(x=x, y=y, pip=1), alpha.realized=0))
#     xy.quant$alpha=alpha
#     xy.quant$xlab=xlab
#     xy.quant$ylab=ylab
#     xy.quant$col=col[1:3]
#     xy.quant$polygon.smooth1=xy.quant$polygon
#     xy.quant$polygon.smooth2=xy.quant$polygon
#     return (xy.quant)
#   }
  
  nknots=round(nknots/2)*2
  
  est.alpha=alpha/nknots
  min.diff=1
  min.diff.est.alpha=est.alpha
  for (iter in 1:40) {
    xy.quant=xy.quantile(x=x, y=y, half.nknots=round(nknots/2), alpha=est.alpha)
#     cat("input=",est.alpha, "  realized=",xy.quant$alpha.realized, "  Target=",alpha, "\n")
    if (abs(alpha-xy.quant$alpha.realized)<=alpha.min.diff) break;
    if (min.diff>=abs(alpha-xy.quant$alpha.realized)) {
      min.diff=abs(alpha-xy.quant$alpha.realized)
      min.diff.est.alpha=est.alpha;
    } else {
      xy.quant=xy.quantile(x=x, y=y, half.nknots=round(nknots/2), alpha=min.diff.est.alpha)
      break
    }
#     est.alpha=min(max(0, est.alpha*alpha/(1e-200+xy.quant$alpha.realized)), 1);
     est.alpha=min(max(0, est.alpha/max(min(xy.quant$alpha.realized/alpha, 100), 0.01)), 1);
    
  }

  xy.quant$alpha=alpha
  xy.quant$xlab=xlab
  xy.quant$ylab=ylab
  xy.quant$col=col[1:3]
  
  class(xy.quant)<-"distfree.cr"
  if (draw) plot(xy.quant)
  xy.quant
}


#obtain angle between (0,0) and (x, y)
xytoangle<-function(x, y) {
  d=sqrt(x^2+y^2)
  ang.cos=acos(x/d)%%(2*pi)
  idx=which(y<0,arr.ind=T)
  ang.cos[idx]=2*pi-ang.cos[idx]
  ang.cos
}

#get index of outer points in x=data.frame(row, col)
border.points.idx<-function(x) {
  x=data.frame(row=x[,1], col=x[,2])

  k=order(x$row)[1]
  bord.points=k;
  ang0=pi;
  
  while(NROW(bord.points)==NROW(unique(bord.points))) { 
    ang=xytoangle(x=x$row-x$row[k], y=x$col-x$col[k])
    len=(x$row-x$row[k])^2+(x$col-x$col[k])^2
    ang_diff=ang-ang0
    idx=which(ang_diff>pi,arr.ind=T) 
    ang_diff=ang_diff%%(2*pi)
    minang=min(ang_diff,  na.rm = T)
    k=which((minang==(ang_diff))&(len>0), arr.ind=T)
    k=k[which(len[k]==max(len[k]))[1]]

    bord.points=c(bord.points, k)
    ang0=ang[k][1]%%(2*pi)
  }
  bord.points
}

