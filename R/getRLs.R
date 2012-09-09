getRLs <-
function(control,treatment){

   #check the feature of raw data (count, log, logratio)
   Check.Log <-function(control){
      logid=c()
      for (i in 1:ncol(control)){
          if (max(control[,i],na.rm=T)<50){
              if (abs(median(control[,i])-1)<0.01)
                  logid[i]='ratio'
              else if (abs(median(control[,i])) < 0.01)
                  logid[i] = 'logratio'
              else 
	          logid[i] = 'log';
           }
           else
               logid[i] = 'count'
       }
       return(logid)
   }
    
   #data preprocessing
   getrank <- function(control,treatment,logid){
      if (logid=='log' | logid=='logratio'){
           control=2^control
           treatment=2^treatment
       }
       control[control<0]=0
       treatment[treatment<0]=0 
       control_treatment=rbind(control,treatment)
       tmp=sort(control_treatment)
       th1=tmp[floor(length(control)/2)]
       th2=th1/10

       control_th=control
       treatment_th=treatment
       control_th[control<th1]=th1
       treatment_th[treatment_th<th1]=th1
       regulate=treatment_th/control_th
       tmp1=sort(regulate,decreasing=T,na.last=F)
       RL=order(regulate,decreasing=T,na.last=F)

       subRL=RL[which(tmp1==1)]
       subcontrol=control[subRL]
       subtreatment=treatment[subRL]
       subcontrol=pmax(subcontrol,th2)
       subtreatment=pmax(subtreatment,th2)
       subregulate=subtreatment/subcontrol
       sub_tmp1=sort(subregulate,decreasing=T)
       indx=order(subregulate,decreasing=T)
       subRL=subRL[indx]

       RL[which(tmp1==1)]=subRL
       RL=order(RL)
    }
    
    probe_num=nrow(control)
    treatment_num=ncol(control)
    RLs=matrix(0,probe_num,treatment_num)
    logid=Check.Log(control)
    for ( i in 1:treatment_num)
        RLs[,i]=getrank(control[,i],treatment[,i],logid[i])
    return(RLs)
}
