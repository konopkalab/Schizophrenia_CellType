# HCP function
hcp <- function(F,Y,k,lambda,lambda2,lambda3,iter=100,tol = 1e-6){
			#### Standardize/Mean center inputs ####
			Yn = scale(Y)
			Fn = scale(F)
			#### Set default values for outputs and initialise B ####
			U = matrix(0,dim(F)[2],k)
			Z = matrix(0,dim(F)[1],k)
			B = matrix(runif(k*dim(Y)[2]),k,dim(Y)[2])
			#### Error checking ####
			if (dim(F)[1] != dim(Y)[1])
  			stop('number of rows in F and Y must agree')
	
			if (k<1 | lambda<tol | lambda2<tol | lambda3<tol)
  			stop('lambda, lambda2, lambda3 must be positive and/or k must be an integer');
			#### Coefficients esitmation ####
			o = matrix(0,1,iter)
			for (ii in 1:iter){
  			o[ii] = norm(Yn-Z %*% B,type='E') + norm(Z-Fn %*% U, type='E')*lambda + norm(B,type='E')*lambda2 + norm(U,type='E')*lambda3
  			Z = ((Yn %*% t(B)) + lambda * (Fn %*% U)) %*% solve((B %*% t(B)) + lambda*diag(dim(B)[1]))
			B = solve(((t(Z) %*% Z) + lambda2 * diag(dim(Z)[2])),(t(Z) %*% Yn))
  			U = solve(((t(Fn) %*% Fn) * lambda + lambda3 * diag(dim(U)[1])),(lambda* t(Fn) %*% Z))
 			if (ii > 1)
    				if ((abs(o[ii]-o[ii-1])/o[ii]) < tol)
				break 
			}
			if (ii>=iter)
  			writeLines("\nPre-mature convergence: Consider increasing the number of iterations")
			#### Error calculation ####
			error = (norm((Yn-Z %*% B),type='E') /norm(Yn,type='E')) + (norm((Z - Fn %*% U),type='E')/norm(Fn%*%U,type='E'))
			error1 = norm(Yn-Z%*%B,type='E')/norm(Yn,type='E')
			error2 = norm(Z-Fn%*%U,type='E')/norm(Fn%*%U,type='E')
			#### Delta change calculation ####
			dz = Z %*% (B %*% t(B) + lambda * diag(dim(B)[1])) - (Yn %*% t(B) + lambda * Fn %*% U)
			db = (t(Z) %*% Z + lambda2 * diag(dim(Z)[2])) %*% B - t(Z) %*% Yn
			du = (t(Fn) %*% Fn *lambda + lambda3 * diag(dim(U)[1])) %*% U - lambda * (t(Fn) %*% Z)
			#### Residual calculation and covariate adjustment ####
			R = (Yn - Z %*% B) * t(replicate(dim(Y)[1],attr(Yn,'scaled:scale'))) + t(replicate(dim(Y)[1],attr(Yn,'scaled:center')))
			#### Final Outputs ####
			return(list(R=R,Z=Z,B=B,U=U,o=o,error=error,error1=error1,error2=error2,dz=dz,db=db,du=du))
}



## VARIANCE EXPLAIN
# counts = gene x sample matrix
# meta = sample x covariate matrix
# threshold = number of PCA to consider (e.g. 5)
# inter = interaction between covariates (e.g. Age:Sex)

VarExp <- function(counts, meta, threshold, inter){
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(optimx))
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)

  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}

  pred.list <- colnames(meta)
  meta <- droplevels(meta)

  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}

  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit,control = lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}

# Bar plot for data visualization
plotVarExp <- function(pvca.res, title){
  suppressPackageStartupMessages(library(ggplot2))
  plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res)
  p <- ggplot2::ggplot(plot.dat, aes(x=eff, y=prop))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::geom_bar(stat="identity", fill="steelblue", colour="steelblue")
  p <- p + ggplot2::geom_text(aes(label=round(prop,3), y=prop+0.04), size=4)
  p <- p + ggplot2::scale_x_discrete(limits=names(pvca.res))
  p <- p + ggplot2::scale_y_continuous(limits = c(0,1))
  p <- p + ggplot2::labs(x= "Effects", y= "WAPV")
  p <- p + ggplot2::theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
}


limmaLOO <- function(mat,pheno,Sample=10,B=5)
                {
                suppressPackageStartupMessages(library(sva))
                suppressPackageStartupMessages(library(limma))
                suppressPackageStartupMessages(library(dplyr))
                suppressPackageStartupMessages(library(tibble))
            
                # Enter the data
                p <- as.data.frame(mat) 
                pd <- pheno
                nLOO=B # Leave multiple out method
                new_pd <- list()
                mod <- list()
                Y.b <- list()
                fitLM <- list()
                fitEb <- list()
                looDGE <- list()
                    for (i in 1:nLOO)
                        {
                            new_pd[[i]] <- pd %>% 
                            rownames_to_column('Sample') %>% 
                            group_by(Diagnosis) %>% 
                            sample_n(Sample) %>% 
                            column_to_rownames('Sample') %>% 
                            as.data.frame() %>%
                            droplevels()
                            Y.b[[i]]=p[,colnames(p) %in% rownames(new_pd[[i]])]
                            Y.b[[i]]=Y.b[[i]][,match(rownames(new_pd[[i]]),colnames(Y.b[[i]]))]
                            mod[[i]] <- model.matrix(~., data = new_pd[[i]])
                            fitLM[[i]] = lmFit(Y.b[[i]],mod[[i]],method="robust");
                            fitEb[[i]] = eBayes(fitLM[[i]]);
                            looDGE[[i]] = topTable(fitEb[[i]], coef = "DiagnosisSchizo",number=nrow(Y.b[[i]]),sort.by="none");
                        }
                        return(looDGE)
                }


limmaPERM <- function(mat,pheno,nPERM=5)
                {
                suppressPackageStartupMessages(library(sva))
                suppressPackageStartupMessages(library(limma))
                suppressPackageStartupMessages(library(parallel))
                p <- as.data.frame(mat) 
                pd <- pheno
                Y.b <- list()
                fitLM <- list()
                fitEb <- list()
                PermDGE <- list()
                    for (i in 1:nPERM)
                        {
                            Y.b <- mclapply(1:nPERM,mc.cores =12, function(i){apply(p, 2, function(col){sample(col)})})
                            colnames(Y.b[[i]]) <- colnames(p)
                            rownames(Y.b[[i]]) <- rownames(p)
                            Y.b[[i]] <- as.data.frame(Y.b[[i]])
                            mod <- model.matrix(~., data = pd)
                            fitLM[[i]] = lmFit(Y.b[[i]],mod,method="robust");
                            fitEb[[i]] = eBayes(fitLM[[i]]);
                            PermDGE[[i]] = topTable(fitEb[[i]], coef = "DiagnosisSchizo",number=nrow(Y.b[[i]]),sort.by="none");
                        }
                        return(PermDGE)
                }





