#PropEdgeUG2D.R;
#Functions for Underlying graphs of PE-PCD in R^2
#################################################################

#' @title The indicator for the presence of an edge from a point to another for the underlying or reflexivity graphs based on
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) - standard basic triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge in the underlying or reflexivity graph of PE-PCDs \eqn{)}
#' for points \code{p1} and \code{p2}.
#'
#' More specifically, it returns
#' \eqn{\max I(}\code{p2} is in \eqn{N_{PE}(p1,r))} \eqn{, I(}\code{p1} is in \eqn{N_{PE}(p2,r))} for points \code{p1} and \code{p2},
#' for the PE-PCD underlying graph.
#' That is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)} or \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise,
#' or
#' \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))} \eqn{\times I(}\code{p1} is in \eqn{N_{PE}(p2,r))} for points \code{p1} and \code{p2},
#' for the PE-PCD reflexivity graph.
#' That is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)} and \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise, where \eqn{N_{PE}(x,r)} is the PE proximity region for point \eqn{x} with expansion parameter \eqn{r \ge 1}.
#'
#' PE proximity region is defined with respect to the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in
#' barycentric coordinates in the interior of the standard basic triangle \eqn{T_b} or based on circumcenter of \eqn{T_b};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
#' \code{rv} is the index of the vertex region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \eqn{T_b}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:stamet2016;textual}{pcds}).
#'
#' @param p1 A 2D point whose PE proximity region is constructed.
#' @param p2 A 2D point. The function determines whether there is an edge from \code{p1} to \code{p1} or not
#' in the underlying and reflexivity graphs of PE-PCDs.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region; must be \eqn{\ge 1}
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle or circumcenter of \eqn{T_b}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
#' @param rv The index of the vertex region in \eqn{T_b} containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#' @param under.graph The type of the graph based on PE-PCDs, "OR" is for the underlying graph, and "AND" is for
#' the reflexivity graph (default is \code{OR}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))} for points \code{p1} and \code{p2}, that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' returns 0 otherwise
#'
#' @seealso \code{\link{IndNPEtri}} and \code{\link{IndNPETe}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C);
#'
#' M<-as.numeric(pcds::runif.bas.tri(1,c1,c2)$g)
#'
#' r<-1.5
#'
#' P1<-as.numeric(pcds::runif.bas.tri(1,c1,c2)$g)
#' P2<-as.numeric(pcds::runif.bas.tri(1,c1,c2)$g)
#' IndUG.PEbas.tri(P1,P2,r,c1,c2,M)
#' IndUG.PEbas.tri(P1,P2,r,c1,c2,M,under.graph = "AND")
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' IndUG.PEbas.tri(P1,P2,r,c1,c2,M)
#' IndUG.PEbas.tri(P2,P1,r,c1,c2,M)
#'
#' #or try
#' Rv<-pcds::rv.bas.tri.cent(P1,c1,c2,M)$rv
#' IndUG.PEbas.tri(P1,P2,r,c1,c2,M,Rv,under.graph = "AND")
#' }
#'
#' @export IndUG.PEbas.tri
IndUG.PEbas.tri <- function(p1,p2,r,c1,c2,M=c(1,1,1),rv=NULL,under.graph="OR")
{
 # if (!is.point(p1) || !is.point(p2))
#  {stop('p1 and p2 must be numeric 2D points')}

 # if (!is.point(r,1) || r<1)
#  {stop('r must be a scalar >= 1')}

 # if (!is.point(c1,1) || !is.point(c2,1))
  #{stop('c1 and c2 must be scalars')}

  #if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  #{stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  #if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  #{stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates
  #        or the circumcenter "CC" ')}

  # y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)
  #
  # CC = pcds::circ.cent.tri(Tb)
  # if (identical(M,"CC") )
  # { M<-CC }
  #
  # if (dimension(M)==3)
  # {M<-bary2cart(M,Tb)}
  #
  # if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,Tb,boundary=FALSE)$in.tri))
  # {stop('center is not the circumcenter or not in the interior of the triangle')}
  #
  # if (isTRUE(all.equal(p1,p2)))
  # {arc<-1; return(arc); stop}
  #
  # if (!in.triangle(p1,Tb,boundary=TRUE)$in.tri || !in.triangle(p2,Tb,boundary=TRUE)$in.tri)
  # {arc<-0; return(arc); stop}
  #
  # if (is.null(rv))
  # { rv<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p1,Tb)$rv,rv.tri.cent(p1,Tb,M)$rv)  #vertex region for p1
  # } else
  # {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  # {stop('vertex index, rv, must be 1, 2 or 3')}}
  #
  # X1<-p1[1]; Y1<-p1[2];
  # X2<-p2[1]; Y2<-p2[2];
  # arc<-0;
  # if (rv==1)
  # {
  #   x1n<-X1*r; y1n<-Y1*r;
  #   if ( Y2 < paraline(c(x1n,y1n),y2,y3,X2)$y ) {arc <-1}
  # } else {
  #   if (rv==2)
  #   {
  #     x1n<-1+(X1-1)*r; y1n<-Y1*r;
  #     if ( Y2 < paraline(c(x1n,y1n),y1,y3,X2)$y ) {arc <-1}
  #   } else {
  #     y1n<-y3[2]+(Y1-y3[2])*r;
  #     if ( Y2 > y1n ) {arc<-1}
  #   }}
  # arc
  arc12 = pcds::IndNPEbas.tri(p1,p2,r,c1,c2,M,rv)
  arc21 = pcds::IndNPEbas.tri(p2,p1,r,c1,c2,M,rv)

  edge = ifelse(under.graph =="OR",max(arc12,arc21),arc12*arc21)
  edge
} #end of the function
#'

