group.block <- function(coords, backg = NULL, k) {
  N <- split(1:nrow(coords), sort(1:nrow(coords)%%k))
  DIST <- spDists(as.matrix(coords)) #create a distance matrix
  rownames(DIST) <- 1:nrow(coords)
  colnames(DIST) <- 1:nrow(coords)
  for(i in 1:k) {
	  if(i == 1) {
      MAX <- which.max(apply(DIST, 1, sum)) #find the furthest point from this group
      DISTmin <- spDistsN1(as.matrix(coords), as.numeric(coords[MAX, ])) #find the distance of all other points from this one
      names(DISTmin) <- rownames(DIST) 
      NAME <- order(DISTmin)[1:length(N[[i]])] 
      GRP <- cbind(coords[NAME, ], group = i)
	    } else {
          DIST <- DIST[-NAME, -NAME]
          coords <- coords[-NAME, ]
          MAX <- which.max(apply(DIST, 1, sum)) 
          DISTmin <- spDistsN1(as.matrix(coords), as.numeric(coords[MAX, ]))
          names(DISTmin) <- rownames(DIST) 
          NAME <- order(DISTmin)[1:length(N[[i]])] 
          GRP <- rbind(GRP, cbind(coords[NAME, ], group = i))
	        }
  }
  GRP <- GRP[order(as.numeric(rownames(GRP))), ]
  if(!is.null(backg)) {
    N <- split(1:nrow(backg), sort(1:nrow(backg)%%k))
    CENT <- foreach(i = 1:k, .combine = "rbind") %do% {
      if(sum(!duplicated(GRP[GRP$group == i, 1:2])) == 1) {
        as.numeric(GRP[GRP$group == i, 1:2][1, ])
      } else {
        as.numeric(centroid(GRP[GRP$group == i, 1:2]))
      }
    }
    DIST <- do.call(cbind.data.frame, lapply(1:k, function(i) spDistsN1(as.matrix(backg), CENT[i, ])))
    POS <- order(apply(DIST, 1, min), decreasing = T)
    GRPB <- NULL
    DIR <- NULL
    for(i in POS) {
      if(!is.null(DIR)) {
        ORD <- order(DIST[i, ])
        ORD <- ORD[-which(ORD %in% DIR)]
        NUM <- ORD[1]
      } else NUM <- order(DIST[i, ])[1]
      GRPB <- rbind(GRPB, cbind(backg[i, ], group = NUM))
      if(nrow(GRPB[GRPB$group == NUM, ]) == length(N[[NUM]])) {
        DIR <- c(DIR, NUM)
      }
    }
    GRPB <- GRPB[order(as.numeric(rownames(GRPB))), ]
  }
  list(GRP, GRPB)
}


