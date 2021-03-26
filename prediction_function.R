Param_opt<-function(dat){

	message("parameter optimizing")
	alpha.runs<-c(0,0.05,0.075,0.1,0.2,0.3,0.4,0.65,0.9,1)

	models.alpha  <- list()
	alpha.vec <- c()

	for (a in 1:length(alpha.runs)){

		message(paste("a=",a,sep=""))

		model=cv.glmnet(as.matrix(dat[,-1]),dat[,cytokine],alpha=alpha.runs[a])

		alpha.vec[a]<- model$cvm[model$lambda == model$lambda.1se]
		models.alpha[[a]] <- model
	}

	message("extract model")

	n=which(alpha.vec == min(alpha.vec))
	best.alpha <- alpha.runs[which(alpha.vec == min(alpha.vec))]
	best.model <- models.alpha[[n]]
	best.lambda<- best.model$lambda.1se

	message(paste("best alpha is",best.alpha))

	message(paste("best lambda is",best.lambda))
}


Model_train<-function(training,testing,best.lambda,best.alpha,result,test){

	modelFit=glmnet(as.matrix(training[,-1]),training[,cytokine],lambda =best.lambda,alpha=best.alpha)
	prediction<-predict(modelFit,newx=as.matrix(testing[,-1]))

	true_v<-testing[,cytokine]
	pred_v<-as.numeric(prediction[,1])

	message(paste("spearman_cor_analysis",cytokine))
	cortest<-cor.test(method="spearman",true_v,pred_v)

	p<-cortest$p.value
        r<-cortest$estimate[[1]]

	result <- list(r=r,p=p)
 	return(result)
}
