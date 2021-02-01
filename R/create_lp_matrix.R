create_lpmatrix <- function(object,newdata,type="lpmatrix",se.fit=FALSE,terms=NULL,exclude=NULL,
                        block.size=NULL,newdata.guaranteed=FALSE,na.action=na.pass,
                        unconditional=FALSE,iterms.type=NULL,...) {
  
 
  if (type!="lpmatrix")  { 
    stop("This function only used for creating lp matrkx")
    type<-"terms"
  }
  
  if (missing(newdata)) na.act <- object$na.action else {
    if (is.null(na.action)) na.act <- NULL 
    else {
      na.txt <- if (is.character(na.action)||is.function(na.action)) get.na.action(na.action) else "na.pass"
      if (na.txt=="na.pass") na.act <- "na.exclude" else
        if (na.txt=="na.exclude") na.act <- "na.omit" else
          na.act <- na.action
    }
  } ## ... done
  
  # get data from which to predict.....  
  nd.is.mf <- FALSE # need to flag if supplied newdata is already a model frame 
  ## get name of response...
  # yname <- all.vars(object$terms)[attr(object$terms,"response")]
  yname <- attr(attr(object$terms,"dataClasses"),"names")[attr(object$terms,"response")]
  if (newdata.guaranteed==FALSE) {
    if (missing(newdata)) { # then "fake" an object suitable for prediction 
      newdata <- object$mf
      new.data.ok <- FALSE
      nd.is.mf <- TRUE
      response <- newdata[[yname]] ## ok even with "cbind(foo,bar)" as yname 
    } else {  # do an R ``standard'' evaluation to pick up data
      new.data.ok <- TRUE
      if (is.data.frame(newdata)&&!is.null(attr(newdata,"terms"))) { # it's a model frame
        if (sum(!(names(object$mf)%in%names(newdata)))) stop(
          "newdata is a model.frame: it should contain all required variables\n")
        nd.is.mf <- TRUE
      } else {
        ## Following is non-standard to allow convenient splitting into blocks
        ## below, and to allow checking that all variables are in newdata ...
        
        ## get names of required variables, less response, but including offset variable
        ## see ?terms.object and ?terms for more information on terms objects
        # yname <- all.vars(object$terms)[attr(object$terms,"response")] ## redundant
        resp <- get.var(yname,newdata,FALSE)
        naresp <- FALSE
        #if (!is.null(object$family$predict)&&!is.null(newdata[[yname]])) {
        if (!is.null(object$family$predict)&&!is.null(resp)) {
          ## response provided, and potentially needed for prediction (e.g. Cox PH family)
          if (!is.null(object$pred.formula)) object$pred.formula <- attr(object$pred.formula,"full")
          response <- TRUE
          Terms <- terms(object)
          #resp <- newdata[[yname]]
          if (is.matrix(resp)) {
            if (sum(is.na(rowSums(resp)))>0) stop("no NAs allowed in response data for this model")
          } else { ## vector response
            if (sum(is.na(resp))>0) {
              naresp <- TRUE ## there are NAs in supplied response
              ## replace them with a numeric code, so that rows are not dropped below
              rar <- range(resp,na.rm=TRUE)
              thresh <- rar[1]*1.01-rar[2]*.01
              resp[is.na(resp)] <- thresh
              newdata[[yname]] <- thresh 
            }
          }  
        } else { ## response not provided
          response <- FALSE 
          Terms <- delete.response(terms(object))
        }
        allNames <- if (is.null(object$pred.formula)) all.vars(Terms) else all.vars(object$pred.formula)
        if (length(allNames) > 0) { 
          ff <- if (is.null(object$pred.formula)) reformulate(allNames) else  object$pred.formula
          if (sum(!(allNames%in%names(newdata)))) { 
            warning("not all required variables have been supplied in  newdata!\n")
          }
          ## note that `xlev' argument not used here, otherwise `as.factor' in 
          ## formula can cause a problem ... levels reset later.
          newdata <- eval(model.frame(ff,data=newdata,na.action=na.act),parent.frame())
          if (naresp) newdata[[yname]][newdata[[yname]]<=thresh] <- NA ## reinstate as NA  
        } ## otherwise it's intercept only and newdata can be left alone
        na.act <- attr(newdata,"na.action")
        #response <- if (response) newdata[[yname]] else NULL
        response <- if (response) get.var(yname,newdata,FALSE) else NULL
      }
    }
  } else { ## newdata.guaranteed == TRUE
    na.act <- NULL
    new.data.ok=TRUE ## it's guaranteed!
    if (!is.null(attr(newdata,"terms"))) nd.is.mf <- TRUE
    #response <- newdata[[yname]]
    response <- get.var(yname,newdata,FALSE)
  }
  
  
  ## now check the factor levels and split into blocks...
  
  if (new.data.ok) {
    ## check factor levels are right ...
    names(newdata)->nn # new data names
    colnames(object$mf)->mn # original names
    for (i in 1:length(newdata)) 
      if (nn[i]%in%mn && is.factor(object$mf[,nn[i]])) { # then so should newdata[[i]] be 
        levm <- levels(object$mf[,nn[i]]) ## original levels
        levn <- levels(factor(newdata[[i]])) ## new levels
        if (sum(!levn%in%levm)>0) { ## check not trying to sneak in new levels 
          msg <- paste("factor levels",paste(levn[!levn%in%levm],collapse=", "),"not in original fit",collapse="")
          warning(msg)
        }
        ## set prediction levels to fit levels...
        if (is.matrix(newdata[[i]])) {
          dum <- factor(newdata[[i]],levels=levm)
          dim(dum) <- dim(newdata[[i]])
          newdata[[i]] <- dum
        } else newdata[[i]] <- factor(newdata[[i]],levels=levm)
      }
    if (type=="newdata") return(newdata)
    
    # split prediction into blocks, to avoid running out of memory
    if (length(newdata)==1) newdata[[2]] <- newdata[[1]] # avoids data frame losing its labels and dimensions below!
    if (is.null(dim(newdata[[1]]))) np <- length(newdata[[1]]) 
    else np <- dim(newdata[[1]])[1] 
    nb <- ncol(object$X)
    if (is.null(block.size)) block.size <- 1000
    if (block.size < 1) block.size <- np
  } else { # no new data, just use object$mf
    np <- nrow(object$mf)
    nb <- ncol(object$X)
  }
  
  if (type=="lpmatrix") block.size <- NULL ## nothing gained by blocking in this case - and offset handling easier this way
  
  ## split prediction into blocks, to avoid running out of memory
  if (is.null(block.size)) { 
    ## use one block as predicting using model frame
    ## and no block size supplied... 
    n.blocks <- 1
    b.size <- array(np,1)
  } else {
    n.blocks <- np %/% block.size
    b.size <- rep(block.size,n.blocks)
    last.block <- np-sum(b.size)
    if (last.block>0) {
      n.blocks <- n.blocks+1  
      b.size[n.blocks] <- last.block
    }
  }
  
  
  # setup prediction arrays...
  ## in multi-linear predictor models, lpi[[i]][j] is the column of model matrix contributing the jth col to lp i 
  lpi <- if (is.list(object$formula)) attr(object$formula,"lpi") else NULL
  nlp <- if (is.null(lpi)) 1 else length(lpi)  ## number of linear predictors
  n.smooth<-length(object$smooth)
  if (type=="lpmatrix") {
    H <- matrix(0,np,nb)
  } else if (type=="terms"||type=="iterms") { 
    term.labels <- attr(object$pterms,"term.labels")
    para.only <- attr(object,"para.only")
    if (is.null(para.only)) para.only <- FALSE  # if TRUE then only return information on parametric part
    n.pterms <- length(term.labels)
    fit <- array(0,c(np,n.pterms+as.numeric(para.only==0)*n.smooth))
    if (se.fit) se <- fit
    ColNames <- term.labels
  } else { ## "response" or "link"
    ## get number of linear predictors, in case it's more than 1...
    #if (is.list(object$formula)) {
    #  nlp <- length(lpi) ## number of linear predictors
    #} else nlp <- 1 
    
    fit <- if (nlp>1) matrix(0,np,nlp) else array(0,np)
    if (se.fit) se <- fit
    fit1 <- NULL ## "response" returned by fam$fv can be non-vector 
  }
  stop <- 0
  if (is.list(object$pterms)) { ## multiple linear predictors
    if (type=="iterms") {
      warning("type iterms not available for multiple predictor cases")
      type <- "terms"
    }
    pstart <- attr(object$nsdf,"pstart") ## starts of parametric blocks in coef vector
    pind <- rep(0,0) ## index of parametric coefs
    Terms <- list();pterms <- object$pterms
    for (i in 1:length(object$nsdf)) {
      Terms[[i]] <- delete.response(object$pterms[[i]])
      if (object$nsdf[i]>0) pind <- c(pind,pstart[i]-1+1:object$nsdf[i])
    }
  } else { ## normal single predictor case
    Terms <- list(delete.response(object$pterms)) ## make into a list anyway
    pterms <- list(object$pterms)
    pstart <- 1
    pind <- 1:object$nsdf ## index of parameteric coefficients
  }
  
  ## check if extended family required intercept to be dropped...
  #drop.intercept <- FALSE 
  #if (!is.null(object$family$drop.intercept)&&object$family$drop.intercept) {
  #  drop.intercept <- TRUE;
  #  ## make sure intercept explicitly included, so it can be cleanly dropped...
  #  for (i in 1:length(Terms)) attr(Terms[[i]],"intercept") <- 1 
  #} 
  drop.intercept <- object$family$drop.intercept
  if (is.null(drop.intercept)) {
    drop.intercept <- rep(FALSE, length(Terms))
  } else {
    ## make sure intercept explicitly included, so it can be cleanly dropped...
    for (i in 1:length(Terms)) {
      if (drop.intercept[i] == TRUE) attr(Terms[[i]],"intercept") <- 1 
    }
  }
  ## index of any parametric terms that have to be dropped
  ## this is used to help with identifiability in multi-
  ## formula models...
  
  drop.ind <- attr(object$nsdf,"drop.ind") 
  
  ####################################
  ## Actual prediction starts here...
  ####################################
  
  s.offset <- NULL # to accumulate any smooth term specific offset
  any.soff <- FALSE # indicator of term specific offset existence
  if (n.blocks > 0) for (b in 1:n.blocks) { # work through prediction blocks
    start <- stop+1
    stop <- start + b.size[b] - 1
    if (n.blocks==1) data <- newdata else data <- newdata[start:stop,]
    X <- matrix(0,b.size[b],nb+length(drop.ind))
    Xoff <- matrix(0,b.size[b],n.smooth) ## term specific offsets 
    offs <- list()
    for (i in 1:length(Terms)) { ## loop for parametric components (1 per lp)
      ## implements safe prediction for parametric part as described in
      ## http://developer.r-project.org/model-fitting-functions.txt
      if (new.data.ok) {
        if (nd.is.mf) mf <- model.frame(data,xlev=object$xlevels) else {
          mf <- model.frame(Terms[[i]],data,xlev=object$xlevels)
          if (!is.null(cl <- attr(pterms[[i]],"dataClasses"))) .checkMFClasses(cl,mf)
        }
        ## next line is just a work around to prevent a spurious warning (e.g. R 3.6) from
        ## model.matrix if contrast relates to a term in mf which is not
        ## part of Terms[[i]] (mode.matrix doc actually defines contrast w.r.t. mf,
        ## not Terms[[i]])...
        oc <- if (length(object$contrasts)==0) object$contrasts else
          object$contrasts[names(object$contrasts)%in%attr(Terms[[i]],"term.labels")]
        Xp <- model.matrix(Terms[[i]],mf,contrasts=oc) 
      } else { 
        Xp <- model.matrix(Terms[[i]],object$mf)
        mf <- newdata # needed in case of offset, below
      }
      offi <- attr(Terms[[i]],"offset")
      if (is.null(offi)) offs[[i]] <- 0 else { ## extract offset
        offs[[i]] <- mf[[names(attr(Terms[[i]],"dataClasses"))[offi+1]]]
      }
      if (drop.intercept[i]) { 
        xat <- attributes(Xp);ind <- xat$assign>0 
        Xp <- Xp[,xat$assign>0,drop=FALSE] ## some extended families need to drop intercept
        xat$assign <- xat$assign[ind];xat$dimnames[[2]]<-xat$dimnames[[2]][ind];
        xat$dim[2] <- xat$dim[2]-1;attributes(Xp) <- xat 
      }
      if (object$nsdf[i]>0) X[,pstart[i]-1 + 1:object$nsdf[i]] <- Xp
    } ## end of parametric loop
    ##   if (length(offs)==1) offs <- offs[[1]] ## messes up later handling
    
    if (!is.null(drop.ind)) X <- X[,-drop.ind]
    
    if (n.smooth) for (k in 1:n.smooth) { ## loop through smooths
      klab <- object$smooth[[k]]$label
      if ((is.null(terms)||(klab%in%terms))&&(is.null(exclude)||!(klab%in%exclude))) {
        Xfrag <- PredictMat(object$smooth[[k]],data)		 
        X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
        Xfrag.off <- attr(Xfrag,"offset") ## any term specific offsets?
        if (!is.null(Xfrag.off)) { Xoff[,k] <- Xfrag.off; any.soff <- TRUE }
      }
      if (type=="terms"||type=="iterms") ColNames[n.pterms+k] <- klab
    } ## smooths done
    
    
    
    if (!is.null(object$Xcentre)) { ## Apply any column centering
      X <- sweep(X,2,object$Xcentre)
    }
    
    # Now have prediction matrix, X, for this block, need to do something with it...
    
    if (type=="lpmatrix") { 
      H[start:stop,] <- X
      if (any.soff) s.offset <- rbind(s.offset,Xoff)
    } else if (type=="terms"||type=="iterms") { ## split results into terms
      lass <- if (is.list(object$assign)) object$assign else list(object$assign)
      k <- 0
      for (j in 1:length(lass)) if (length(lass[[j]])) { ## work through assign list
        ind <- 1:length(lass[[j]]) ## index vector for coefs involved
        nptj <- max(lass[[j]]) ## number of terms involved here
        if (nptj>0) for (i in 1:nptj) { ## work through parametric part
          k <- k + 1 ## counts total number of parametric terms
          ii <- ind[lass[[j]]==i] + pstart[j] - 1 
          fit[start:stop,k] <- X[,ii,drop=FALSE]%*%object$coefficients[ii]
          if (se.fit) se[start:stop,k] <-
            sqrt(pmax(0,rowSums((X[,ii,drop=FALSE]%*%object$Vp[ii,ii])*X[,ii,drop=FALSE])))
        }
      } ## assign list done
      if (n.smooth&&!para.only) {
        for (k in 1:n.smooth) # work through the smooth terms 
        { first <- object$smooth[[k]]$first.para; last <- object$smooth[[k]]$last.para
        fit[start:stop,n.pterms+k] <- X[,first:last,drop=FALSE] %*% object$coefficients[first:last] + Xoff[,k]
        if (se.fit) { # diag(Z%*%V%*%t(Z))^0.5; Z=X[,first:last]; V is sub-matrix of Vp
          if (type=="iterms"&& attr(object$smooth[[k]],"nCons")>0) { ## termwise se to "carry the intercept
            ## some general families, add parameters after cmX created, which are irrelevant to cmX... 
            if (length(object$cmX) < ncol(X)) object$cmX <- c(object$cmX,rep(0,ncol(X)-length(object$cmX)))
            if (!is.null(iterms.type)&&iterms.type==2) object$cmX[-(1:object$nsdf)] <- 0 ## variability of fixed effects mean only
            X1 <- matrix(object$cmX,nrow(X),ncol(X),byrow=TRUE)
            meanL1 <- object$smooth[[k]]$meanL1
            if (!is.null(meanL1)) X1 <- X1 / meanL1              
            X1[,first:last] <- X[,first:last]
            se[start:stop,n.pterms+k] <- sqrt(pmax(0,rowSums((X1%*%object$Vp)*X1)))
          } else se[start:stop,n.pterms+k] <- ## terms strictly centred
              sqrt(pmax(0,rowSums((X[,first:last,drop=FALSE]%*%
                                     object$Vp[first:last,first:last,drop=FALSE])*X[,first:last,drop=FALSE])))
        } ## end if (se.fit)
        } 
        colnames(fit) <- ColNames
        if (se.fit) colnames(se) <- ColNames
      } else {
        if (para.only&&is.list(object$pterms)) { 
          ## have to use term labels that match original data, or termplot fails 
          ## to plot. This only applies for 'para.only==1' calls which are 
          ## designed for use from termplot called from plot.gam
          term.labels <- unlist(lapply(object$pterms,attr,"term.labels"))
        }
        colnames(fit) <- term.labels
        if (se.fit) colnames(se) <- term.labels
        if (para.only) { 
          # retain only terms of order 1 - this is to make termplot work
          order <- if (is.list(object$pterms)) unlist(lapply(object$pterms,attr,"order")) else attr(object$pterms,"order")
          term.labels <- term.labels[order==1]
          ## fit <- as.matrix(as.matrix(fit)[,order==1])
          fit <- fit[,order==1,drop=FALSE]
          colnames(fit) <- term.labels
          if (se.fit) { ## se <- as.matrix(as.matrix(se)[,order==1])
            se <- se[,order==1,drop=FALSE]
            colnames(se) <- term.labels 
          }
        } 
      } 
    } else { ## "link" or "response" case
      fam <- object$family
      k <- attr(attr(object$mf,"terms"),"offset")
      if (nlp>1) { ## multiple linear predictor case
        if (is.null(fam$predict)||type=="link") {
          ##pstart <- c(pstart,ncol(X)+1)
          ## get index of smooths with an offset...
          off.ind <- (1:n.smooth)[as.logical(colSums(abs(Xoff)))]
          for (j in 1:nlp) { ## looping over the model formulae
            ind <- lpi[[j]] ##pstart[j]:(pstart[j+1]-1)
            fit[start:stop,j] <- X[,ind,drop=FALSE]%*%object$coefficients[ind] + offs[[j]]
            if (length(off.ind)) for (i in off.ind) { ## add any term specific offsets
              if (object$smooth[[i]]$first.para%in%ind)  fit[start:stop,j] <- fit[start:stop,j] + Xoff[,i]
            }
            if (se.fit) se[start:stop,j] <- 
              sqrt(pmax(0,rowSums((X[,ind,drop=FALSE]%*%object$Vp[ind,ind,drop=FALSE])*X[,ind,drop=FALSE])))
            ## model offset only handled for first predictor... fixed
            ##if (j==1&&!is.null(k))  fit[start:stop,j] <- fit[start:stop,j] + model.offset(mf)
            if (type=="response") { ## need to transform lp to response scale
              linfo <- object$family$linfo[[j]] ## link information
              if (se.fit) se[start:stop,j] <- se[start:stop,j]*abs(linfo$mu.eta(fit[start:stop,j]))
              fit[start:stop,j] <- linfo$linkinv(fit[start:stop,j])
            }
          } ## end of lp loop
        } else { ## response case with own predict code
          #lpi <- list();pst <- c(pstart,ncol(X)+1)
          #for (i in 1:(length(pst)-1)) lpi[[i]] <- pst[i]:(pst[i+1]-1)
          attr(X,"lpi") <- lpi  
          ffv <- fam$predict(fam,se.fit,y=response[start:stop],X=X,beta=object$coefficients,
                             off=offs,Vb=object$Vp)
          if (is.matrix(fit)&&!is.matrix(ffv[[1]])) {
            fit <- fit[,1]; if (se.fit) se <- se[,1]
          }
          if (is.matrix(ffv[[1]])&&(!is.matrix(fit)||ncol(ffv[[1]])!=ncol(fit))) {
            fit <- matrix(0,np,ncol(ffv[[1]])); if (se.fit) se <- fit
          }
          if (is.matrix(fit)) {
            fit[start:stop,] <- ffv[[1]]
            if (se.fit) se[start:stop,] <- ffv[[2]]
          } else {
            fit[start:stop] <- ffv[[1]]
            if (se.fit) se[start:stop] <- ffv[[2]]
          }
        } ## end of own response prediction code
      } else { ## single linear predictor
        offs <- if (is.null(k)) rowSums(Xoff) else rowSums(Xoff) + model.offset(mf)
        fit[start:stop] <- X%*%object$coefficients + offs
        if (se.fit) se[start:stop] <- sqrt(pmax(0,rowSums((X%*%object$Vp)*X)))
        if (type=="response") { # transform    
          linkinv <- fam$linkinv
          if (is.null(fam$predict)) {
            dmu.deta <- fam$mu.eta  
            if (se.fit) se[start:stop]<-se[start:stop]*abs(dmu.deta(fit[start:stop])) 
            fit[start:stop] <- linkinv(fit[start:stop])
          } else { ## family has its own prediction code for response case
            ffv <- fam$predict(fam,se.fit,y=response[start:stop],X=X,beta=object$coefficients,off=offs,Vb=object$Vp)
            if (is.null(fit1)&&is.matrix(ffv[[1]])) {
              fit1 <- matrix(0,np,ncol(ffv[[1]]))
              if (se.fit) se1 <- fit1
            }
            if (is.null(fit1)) {
              fit[start:stop] <- ffv[[1]]
              if (se.fit) se[start:stop] <- ffv[[2]]
            } else {
              fit1[start:stop,] <- ffv[[1]]
              if (se.fit) se1[start:stop,] <- ffv[[2]]
            }
          }
        }
      } ## single lp done
    } ## end of link or response case 
    rm(X)
  } ## end of prediction block loop
  
  if ((type=="terms"||type=="iterms")&&(!is.null(terms)||!is.null(exclude))) { # return only terms requested via `terms'
    cnames <- colnames(fit)
    if (!is.null(terms)) {
      if (sum(!(terms %in%cnames))) 
        warning("non-existent terms requested - ignoring")
      else { 
        fit <- fit[,terms,drop=FALSE]
        if (se.fit) {
          se <- se[,terms,drop=FALSE]
        }
      }
    }
    if (!is.null(exclude)) {
      if (sum(!(exclude %in%cnames))) 
        warning("non-existent exclude terms requested - ignoring")
      else { 
        exclude <- which(cnames%in%exclude) ## convert to numeric column index
        fit <- fit[,-exclude,drop=FALSE]
        if (se.fit) {
          se <- se[,-exclude,drop=FALSE]
        }
      }
    }
  }
  
  if (type=="response"&&!is.null(fit1)) {
    fit <- fit1
    if (se.fit) se <- se1
  }
  
  rn <- rownames(newdata)
  if (type=="lpmatrix") { 
    # TODO: need to find a way to make the names
    colnames(H) <- G$term.names;rownames(H)<-rn
    if (!is.null(s.offset)) { 
      s.offset <- napredict(na.act,s.offset)
      attr(H,"offset") <- s.offset ## term specific offsets...
    }
    #if (!is.null(attr(attr(object$mf,"terms"),"offset"))) {
    #  attr(H,"model.offset") <- napredict(na.act,model.offset(mf)) 
    #}
    if (!is.null(offs)) {
      offs <- offs[1:nlp]
      for (i in 1:nlp) offs[[i]] <- napredict(na.act,offs[[i]])
      attr(H,"model.offset") <- if (nlp==1) offs[[1]] else offs
    }
    H <- napredict(na.act,H)
    if (length(object$nsdf)>1) { ## add "lpi" attribute if more than one l.p.
      #lpi <- list();pst <- c(pstart,ncol(H)+1)
      #for (i in 1:(length(pst)-1)) lpi[[i]] <- pst[i]:(pst[i+1]-1)  
      attr(H,"lpi") <- lpi
    }
  } else { 
    if (se.fit) { 
      if (is.null(nrow(fit))) {
        names(fit) <- rn
        names(se) <- rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se) 
      } else { 
        rownames(fit)<-rn
        rownames(se)<-rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se)
      }
      H<-list(fit=fit,se.fit=se) 
    } else { 
      H <- fit
      if (is.null(nrow(H))) names(H) <- rn else
        rownames(H)<-rn
      H <- napredict(na.act,H)
    }
  }
  if ((type=="terms"||type=="iterms")&&attr(object$terms,"intercept")==1) attr(H,"constant") <- object$coefficients[1]
  H # ... and return
} ## end of predict.gam

