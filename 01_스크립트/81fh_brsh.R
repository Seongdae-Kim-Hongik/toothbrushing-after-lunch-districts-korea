# 사용: ST_GRAPH=<593 그래프> ST_THREADS=1 rslot.sh Rscript 81fh_fit_fayherriot.R <TAG>   (TAG = DXAW_H | DXAW_D, ST_VER=std3 패널 사용)
#  주모형: y=logit(SR/100) ~ N(η, v) [v=SE²/(p(1-p))² 고정, gaussian scale=1/v, 관측 정밀도 초모수 고정]
#          η = 1 + forced + 가설 + f(area,bym2) + f(time,rw1) + f(yiid,iid) + f(id,iid)[FH 모형 오차]
#  단계(각 CSV, 있으면 건너뜀): S0 공간 선별 지표 · S2 세트 비교 · S3 주모형·음성대조·공간 OR·잔차 Moran·분산 몫·강제통제만 대비
#   · S4 창별+상호작용 · S5 COVID · S6 연도별 Moran · S7 다중우주 · S8 KH II/III/IV · S9 지역 간·내 분해 · S10 지연 1·2년
#   · S11 시도 하나씩 빼기 · S12 우도 4종(이항 n_eff·베타·무가중 가우시안·조율 CR+65+) · S13 queen 그래프 · S14 격차 G(2차 결과) · S15 PC prior 대안
#  설계 근거: 04_분석/DXAW/00_DXAW_설계메모_적합전고정_261110.md
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);TAG<-args[1];stopifnot(TAG=="BRSH")
ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
SHP<-Sys.getenv("ST_SHP",file.path(normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),"..")),"data","boundary","final.shp"))
SUF<-"_std3";OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,"results_std3");dir.create(RES,showWarnings=FALSE)
NT<-as.integer(Sys.getenv("ST_THREADS",unset="1"));inla.setOption(num.threads=sprintf("%d:1",NT))
ONLY<-Sys.getenv("ST_FH_ONLY",unset="");RUN<-function(s)!nzchar(ONLY)||s%in%strsplit(ONLY,",")[[1]]
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE);RU<-S$rules
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
META<-read.csv(file.path(OUTD,sprintf("panel_META_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
LAB<-setNames(META$label_en,META$variable);LAB["ext__dental_care_sir"]<-"Dental care use, age-standardised SIR (log)";LAB["dentists_per_100k"]<-"Dentists per 100,000 (log)";LAB["ext__manufacturing_share"]<-"Manufacturing employment share (%)"
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
FO<-read.csv(file.path(ROOT,"02_패널","BRSH_fh_outcomes.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);FO$fbd_region<-enc2utf8(FO$fbd_region)
DD<-merge(DD,FO[,c("year","fbd_region","SR","SR_SE","N","CR","CR_SE","y_logit","v_logit","y_logit_cr","v_logit_cr","n_eff","y_eff")],by=c("year","fbd_region"),all.x=TRUE,sort=FALSE)
stopifnot(nrow(DD)==229*11,all(!is.na(DD$y_logit)),max(abs(DD$SR-DD$cases))<1e-9)
FORCED<-S$forced;NEG<-S$negative_controls;SETS<-lapply(S$sets_after_vif,as.character);HYPV<-SETS[[RU$main_set]];MAINV<-unique(c(HYPV,FORCED))
CAND<-setdiff(unique(unlist(SETS[setdiff(names(SETS),"S1b_obesity_sensitivity")])),FORCED);stopifnot(!("obesity_rate"%in%CAND));NMV<-RU$multiverse_n;if(nzchar(Sys.getenv("ST_MV_N")))NMV<-as.integer(Sys.getenv("ST_MV_N"))
WIN<-INFO$windows;POSTY<-RU$covid_post_from
LOG<-function(...)cat(sprintf("[81fh %s %s] ",TAG,format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
stopifnot(all(c(MAINV,NEG)%in%names(DD)))

## ── 그래프(593 외부 그래프, 81과 같은 방식; 그래프 파일은 tempdir·pid별 — 병렬 작업 동시 쓰기 경합 방지 11/10) ──
G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);stopifnot(all(unique(DD$fbd_region)%in%G$region))
nbq<-poly2nb(G,queen=TRUE);isl<-which(card(nbq)==0)
if(length(isl)){cen<-st_coordinates(st_centroid(st_geometry(G)));kn<-knearneigh(cen,k=1)$nn
  for(i in isl){j<-kn[i,1];nbq[[i]]<-sort(unique(c(nbq[[i]][nbq[[i]]>0],j)));nbq[[j]]<-sort(unique(c(nbq[[j]][nbq[[j]]>0],i)))}}
nb<-nbq;GP<-Sys.getenv("ST_GRAPH");stopifnot(nzchar(GP))
sg<-inla.read.graph(GP);stopifnot(sg$n==nrow(G));for(i in seq_len(sg$n))nb[[i]]<-as.integer(sort(sg$nbs[[i]]));attr(nb,"sym")<-TRUE
GF<-file.path(tempdir(),sprintf("nb_%s_%d.graph",TAG,Sys.getpid()));nb2INLA(GF,nb);g<-inla.read.graph(GF);lw<-nb2listw(nb,style="W",zero.policy=TRUE)
GFQ<-file.path(tempdir(),sprintf("nbqueen_%s_%d.graph",TAG,Sys.getpid()));nb2INLA(GFQ,nbq);gq<-inla.read.graph(GFQ)
n<-nrow(G);mI<-function(x,alt="greater"){mm<-moran.mc(x,lw,999,zero.policy=TRUE,alternative=alt);c(I=unname(mm$statistic),p=mm$p.value,EI=-1/(length(x)-1))}   # 2026-11-18 잔차는 alt='two.sided'

## ── 색인 ──
DD$area<-match(DD$fbd_region,G$region);DD<-DD[order(DD$year,DD$area),];Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$id<-seq_len(nrow(DD))
DD$w<-1/DD$v_logit;DD$w_cr<-1/DD$v_logit_cr;DD$p<-DD$SR/100
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
bt<-function(v)sprintf("`%s`",v)
SPf<-function(gr="g",hp="pc")sprintf("f(area,model='bym2',graph=%s,scale.model=TRUE,constr=TRUE,hyper=%s)",gr,hp)
TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)";CELL<-"f(id,model='iid',hyper=pcp)"
FIXP<-list(hyper=list(prec=list(initial=0,fixed=TRUE)))
fitfh<-function(vars,data=DD,lik="fh",extra=character(0),re=TRUE,cell=TRUE,cpo=TRUE,sp=SPf(),robust=FALSE){
  yv<-switch(lik,fh="y_logit",gfree="y_logit",cr="y_logit_cr",meas="y_logit_meas",binom="y_eff",beta="p")
  usecell<-cell&&lik%in%c("fh","cr","meas","binom")
  TIv<-sprintf("f(time,model='rw1',scale.model=TRUE,hyper=pcp,values=1:%d)",max(data$time))
  rhs<-c("1",bt(vars),extra,if(re)c(sp,TIv,IT),if(usecell)CELL);form<-as.formula(paste(yv,"~",paste(rhs,collapse="+")))
  ci<-if(robust)list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian") else list(control.vb=list(emergency=50))
  cc<-list(dic=TRUE,waic=TRUE,cpo=cpo,config=FALSE)
  if(lik%in%c("fh","cr","meas")){sc<-switch(lik,fh=data$w,cr=data$w_cr,meas=1/data$v_logit_meas)
    return(inla(form,family="gaussian",data=data,scale=sc,control.family=FIXP,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))}
  if(lik=="gfree")return(inla(form,family="gaussian",data=data,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))
  if(lik=="binom")return(inla(form,family="binomial",Ntrials=data$n_eff,data=data,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))
  if(lik=="beta")return(inla(form,family="beta",data=data,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))
  stop(lik)}
safe<-function(...)tryCatch(fitfh(...),error=function(e){message("  ✗ fit failed: ",conditionMessage(e)," → retry robust");tryCatch(fitfh(...,robust=TRUE),error=function(e2){message("  ✗ retry failed: ",conditionMessage(e2));NULL})})
orr<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];data.frame(cov=v,OR=round(exp(r$mean),4),lo=round(exp(r$`0.025quant`),4),hi=round(exp(r$`0.975quant`),4),
  P_gt1=round(1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]]),4),cred=(r$`0.025quant`>0|r$`0.975quant`<0),stringsAsFactors=FALSE)}
fstat<-function(m){cp<-m$cpo$cpo;hp<-m$summary.hyperpar;data.frame(WAIC=round(m$waic$waic,1),DIC=round(m$dic$dic,1),LCPO=if(is.null(cp))NA else round(-mean(log(cp[is.finite(cp)&cp>0])),4),
  phi=round(if("Phi for area"%in%rownames(hp))hp["Phi for area","mean"] else NA,3),
  sd_area=round(if("Precision for area"%in%rownames(hp))1/sqrt(hp["Precision for area","mean"]) else NA,4),sd_cell=round(if("Precision for id"%in%rownames(hp))1/sqrt(hp["Precision for id","mean"]) else NA,4))}
tab<-function(m,vars,extra=list()){r<-do.call(rbind,lapply(vars,function(v)orr(m,v)));r$label_en<-LAB[r$cov];r$role<-ifelse(r$cov%in%FORCED,"forced","hypothesis");for(k in names(extra))r[[k]]<-extra[[k]];r}
timed<-function(label,expr){t0<-Sys.time();r<-expr;LOG("%s · %.1f min",label,as.numeric(difftime(Sys.time(),t0,units="mins")));r}
clusters<-function(nodes){nodes<-as.integer(nodes);seen<-integer(0);best<-0;for(s in nodes){if(s%in%seen)next;st<-s;comp<-integer(0)
  while(length(st)){v<-st[1];st<-st[-1];if(v%in%seen)next;seen<-c(seen,v);comp<-c(comp,v);st<-c(st,setdiff(intersect(nb[[v]],nodes),seen))};best<-max(best,length(comp))};best}
sidoOf<-sapply(G$region,function(r)sub("^(서울시|부산시|대구시|인천시|광주시|대전시|울산시|세종시|경기도|강원도|충청북도|충청남도|전라북도|전라남도|경상북도|경상남도|제주도).*$","\\1",r))
spatialOR<-function(m){ur<-m$summary.random$area[1:n,];pep<-sapply(m$marginals.random$area[1:n],function(mg)1-inla.pmarginal(0,mg))
  rr<-data.frame(area=1:n,fbd_region=G$region,sido=sidoOf,OR=round(exp(ur$mean),4),lo=round(exp(ur$`0.025quant`),4),hi=round(exp(ur$`0.975quant`),4),P_gt1=round(pep,4))
  rr$class<-ifelse(rr$P_gt1>0.95,"hotspot",ifelse(rr$P_gt1<0.05,"coldspot","neutral"));rr}

## ── S0 공간 선별 지표(설계 메모 §4) ──
F0<-file.path(RES,"S0_screen.csv")
if(RUN("S0")&&!file.exists(F0)){set.seed(20261110)
  pooled<-tapply(DD$y_logit,DD$area,mean);mo<-mI(as.numeric(pooled))
  yI<-sapply(sort(unique(DD$year)),function(y){d<-DD[DD$year==y,];x<-d$y_logit[order(d$area)];mI(x)["I"]})
  srA<-tapply(DD$SR[DD$year<=2016],DD$area[DD$year<=2016],mean);srB<-tapply(DD$SR[DD$year>=2022],DD$area[DD$year>=2022],mean);pers<-cor(srA,srB,method="spearman")
  srP<-tapply(DD$SR,DD$area,mean);p9010<-as.numeric(quantile(srP,.9)/quantile(srP,.1))
  m<-timed("S0 forced-only",safe(FORCED));fs<-fstat(m);rr<-spatialOR(m)
  out<-data.frame(tag=TAG,moran_I_pooled_logitSR=round(mo["I"],4),moran_p=mo["p"],moran_I_year_min=round(min(yI),4),moran_I_year_max=round(max(yI),4),persistence_rank_rho=round(pers,3),
    P90_P10_SR=round(p9010,3),phi_forced_only=fs$phi,n_hot=sum(rr$class=="hotspot"),n_cold=sum(rr$class=="coldspot"),largest_hot_cluster=clusters(rr$area[rr$class=="hotspot"]),largest_cold_cluster=clusters(rr$area[rr$class=="coldspot"]))
  out$track<-ifelse(out$moran_p<0.01,ifelse(out$moran_I_pooled_logitSR<0.45&out$moran_I_pooled_logitSR>=0.30,"공간(약한 공간형)",ifelse(out$moran_I_pooled_logitSR>=0.45,"공간","공간(Moran<0.30, 약신호)")),"도농 트랙(Moran 비유의)")
  write.csv(out,F0,row.names=FALSE);write.csv(rr,file.path(RES,"S0_forced_only_spatialOR.csv"),row.names=FALSE);LOG("S0: Moran %.3f (p %.3f) phi %.2f track %s",out$moran_I_pooled_logitSR,out$moran_p,out$phi_forced_only,out$track);rm(m)}

## ── S2 세트 비교 ──
F2<-file.path(RES,"S2_setcompare.csv");F2e<-file.path(RES,"S2_setcompare_OR.csv")
if(RUN("S2")&&!file.exists(F2)){cmp<-NULL;est<-NULL
  for(sn in c("REF_forced_only",names(SETS))){cv<-if(sn=="REF_forced_only")FORCED else unique(c(SETS[[sn]],FORCED));m<-timed(sprintf("S2 %s",sn),safe(cv));if(is.null(m))next
    cmp<-rbind(cmp,cbind(set=sn,nvar=length(cv),n_hyp=length(cv)-length(FORCED),fstat(m)));if(sn!="REF_forced_only")est<-rbind(est,tab(m,cv,list(set=sn)))}
  write.csv(cmp,F2,row.names=FALSE);write.csv(est,F2e,row.names=FALSE)}

## ── S3 주모형 ──
F3<-file.path(RES,"S3_main_OR.csv");F3f<-file.path(RES,"S3_main_fit.csv");F3n<-file.path(RES,"S3_negcontrol.csv");F3s<-file.path(RES,"S3_spatial_OR.csv");F3v<-file.path(RES,"S3_varpart.csv");F3y<-file.path(RES,"S3_year_effect.csv")
if(RUN("S3")&&!all(file.exists(F3,F3f,F3n,F3s,F3v))){
  m<-timed("S3 main",safe(MAINV));stopifnot(!is.null(m));ir<-tab(m,MAINV);ir$expected<-sapply(ir$cov,function(v){x<-S$vars[[v]]$expected_direction;if(is.null(x))"" else substr(x,1,1)})
  write.csv(ir,F3,row.names=FALSE)
  hp<-m$summary.hyperpar;write.csv(data.frame(hyperpar=rownames(hp),round(hp[,c("mean","0.025quant","0.975quant")],4)),file.path(RES,"S3_main_hyper.csv"),row.names=FALSE)
  rr<-spatialOR(m);agg<-aggregate(cbind(SR,y_logit)~area,DD,mean);rr$SR_mean<-round(agg$SR[match(rr$area,agg$area)],2);write.csv(rr,F3s,row.names=FALSE)
  ye<-m$summary.random$time;write.csv(data.frame(year=Y0+ye$ID-1,rw1=round(ye$mean,4),lo=round(ye$`0.025quant`,4),hi=round(ye$`0.975quant`,4)),F3y,row.names=FALSE)
  # 잔차 = y − (선형예측 − 셀 iid) 를 시군구 평균 → Moran (FH 모형 오차 포함 잔차)
  eta<-m$summary.fitted.values$mean;cellm<-m$summary.random$id$mean[DD$id];res<-DD$y_logit-(eta-cellm);ra<-tapply(res,DD$area,mean);set.seed(3);mr<-mI(as.numeric(ra),alt="two.sided")
  mo<-mI(as.numeric(tapply(DD$y_logit,DD$area,mean)))
  # 분산 몫(사후평균 성분의 분산): 고정효과 X·β, 공간, 시간(rw1+iid), 셀
  X<-as.matrix(DD[,MAINV]);bfix<-m$summary.fixed[MAINV,"mean"];xb<-X%*%bfix;xbH<-as.matrix(DD[,HYPV])%*%m$summary.fixed[HYPV,"mean"];xbF<-as.matrix(DD[,FORCED])%*%m$summary.fixed[FORCED,"mean"]
  spc<-m$summary.random$area$mean[DD$area];tim<-m$summary.random$time$mean[match(DD$time,m$summary.random$time$ID)]+m$summary.random$yiid$mean[match(DD$yiid,m$summary.random$yiid$ID)]
  vv<-c(hypothesis=var(xbH),forced=var(xbF),spatial=var(spc),temporal=var(tim),cell=var(cellm),sampling=mean(DD$v_logit));write.csv(data.frame(component=names(vv),variance=round(vv,6),share=round(vv/sum(vv),4)),F3v,row.names=FALSE)
  mf<-timed("S3 forced-only for comparison",safe(FORCED));rf<-spatialOR(mf)
  cmpF<-data.frame(spearman_spatial_OR=round(cor(log(rr$OR),log(rf$OR),method="spearman"),4),sd_logOR_main=round(sd(log(rr$OR)),4),sd_logOR_forced=round(sd(log(rf$OR)),4),
    hot_main=sum(rr$class=="hotspot"),hot_forced=sum(rf$class=="hotspot"),cold_main=sum(rr$class=="coldspot"),cold_forced=sum(rf$class=="coldspot"),WAIC_main=round(m$waic$waic,1),WAIC_forced=round(mf$waic$waic,1))
  write.csv(cmpF,file.path(RES,"S3_vs_forced_only.csv"),row.names=FALSE)
  fs<-fstat(m);write.csv(data.frame(fs,n_obs=nrow(DD),n_regions=n,n_years=11,n_hotspot=sum(rr$class=="hotspot"),n_coldspot=sum(rr$class=="coldspot"),OR_min=min(rr$OR),OR_max=max(rr$OR),
    n_cred_hyp=sum(ir$cred[ir$role=="hypothesis"]),cred_hyp=paste(ir$cov[ir$cred&ir$role=="hypothesis"],collapse=";"),obs_moran_I=round(mo["I"],4),obs_moran_p=mo["p"],resid_moran_I=round(mr["I"],4),resid_moran_p=mr["p"],
    largest_hot=clusters(rr$area[rr$class=="hotspot"]),largest_cold=clusters(rr$area[rr$class=="coldspot"])),F3f,row.names=FALSE)
  nc<-NULL;for(v in NEG){mm<-timed(sprintf("S3 negctrl %s",v),safe(c(MAINV,v)));if(is.null(mm))next;r<-orr(mm,v);r$label_en<-LAB[v];r$WAIC<-round(mm$waic$waic,1);nc<-rbind(nc,r)};write.csv(nc,F3n,row.names=FALSE)
  LOG("S3: cred hyp = %s; hot %d cold %d; resid Moran %.3f",paste(ir$cov[ir$cred&ir$role=="hypothesis"],collapse=","),sum(rr$class=="hotspot"),sum(rr$class=="coldspot"),mr["I"]);rm(m,mf)}

## ── S4 창별 + 상호작용 ──
F4<-file.path(RES,"S4_window_OR.csv");F4i<-file.path(RES,"S4_interaction.csv")
if(RUN("S4")&&!all(file.exists(F4,F4i))){wr<-NULL
  for(wn in names(WIN)){yy<-WIN[[wn]];d<-DD[DD$year>=yy[1]&DD$year<=yy[2],];d$time<-d$year-min(d$year)+1;d$yiid<-d$time;d$id<-seq_len(nrow(d))
    m<-timed(sprintf("S4 %s",wn),safe(MAINV,data=d));if(is.null(m))next;wr<-rbind(wr,tab(m,MAINV,list(window=wn,WAIC=round(m$waic$waic,1))))}
  write.csv(wr,F4,row.names=FALSE)
  d<-DD;wl<-names(WIN);d$W<-factor(d$window,levels=wl);ex<-character(0)
  for(k in 2:length(wl)){dn<-sprintf("W%d",k);d[[dn]]<-as.integer(d$W==wl[k]);ex<-c(ex,dn);for(v in HYPV){nm<-sprintf("%s_x_W%d",v,k);d[[nm]]<-d[[v]]*d[[dn]];ex<-c(ex,bt(nm))}}
  m<-timed("S4 interaction",safe(MAINV,data=d,extra=ex));ir<-NULL;for(k in 2:length(wl))for(v in HYPV){r<-orr(m,sprintf("%s_x_W%d",v,k));r$cov<-v;r$window<-wl[k];ir<-rbind(ir,r)};write.csv(ir,F4i,row.names=FALSE);rm(m)}

## ── S5 COVID ──
F5<-file.path(RES,"S5_covid.csv")
if(RUN("S5")&&!file.exists(F5)){d<-DD;ex<-"post";for(v in HYPV){nm<-sprintf("%s_x_post",v);d[[nm]]<-d[[v]]*d$post;ex<-c(ex,bt(nm))}
  m<-timed("S5 covid",safe(MAINV,data=d,extra=ex));out<-NULL;for(v in HYPV){r<-orr(m,sprintf("%s_x_post",v));r$cov<-v;r$term<-"var_x_post";out<-rbind(out,r)};write.csv(out,F5,row.names=FALSE);rm(m)}

## ── S6 연도별 Moran ──
F6<-file.path(RES,"S6_moran_year.csv")
if(RUN("S6")&&!file.exists(F6)){set.seed(2);my<-do.call(rbind,lapply(sort(unique(DD$year)),function(y){d<-DD[DD$year==y,];x<-d$y_logit[order(d$area)];mm<-mI(x)
  data.frame(year=y,moran_I=round(mm["I"],4),p=mm["p"],SR_popweighted=round(weighted.mean(d$SR,d$population),2))}));write.csv(my,F6,row.names=FALSE)
  ken<-NA;try({R<-sapply(sort(unique(DD$year)),function(y)rank(DD$SR[DD$year==y][order(DD$area[DD$year==y])]));m_<-ncol(R);nn<-nrow(R);Rs<-rowSums(R);ken<-12*sum((Rs-mean(Rs))^2)/(m_^2*(nn^3-nn))})
  write.csv(data.frame(kendall_W=round(ken,3)),file.path(RES,"S6_kendallW.csv"),row.names=FALSE)}

## ── S7 다중우주 ──
F7s<-file.path(RES,"S7_multiverse_specs.csv");F7r<-file.path(RES,"S7_multiverse_results.csv");F7<-file.path(RES,"S7_multiverse_summary.csv")
if(RUN("S7")&&!file.exists(F7)){set.seed(20261110);POOL<-CAND
  if(!file.exists(F7s)){specs<-list();keys<-character(0);while(length(specs)<NMV){inc<-POOL[runif(length(POOL))<0.5];if(length(inc)<1)next;k<-paste(sort(inc),collapse="|");if(k%in%keys)next;specs[[length(specs)+1]]<-inc;keys<-c(keys,k)}
    write.csv(data.frame(spec=seq_along(specs),nvar=sapply(specs,length),vars=sapply(specs,paste,collapse="|")),F7s,row.names=FALSE)}
  if(Sys.getenv("ST_MV_STOP")=="1"){LOG("S7 specs written");quit(save="no")}
  RNG<-Sys.getenv("ST_MV_RANGE");SPSEL<-NULL;if(nzchar(RNG)){ab<-as.integer(strsplit(RNG,"-")[[1]]);SPSEL<-ab[1]:ab[2];F7r<-file.path(RES,sprintf("S7_multiverse_results_part_%d_%d.csv",ab[1],ab[2]))}
  PARTS<-list.files(RES,pattern="^S7_multiverse_results_part_.*\\.csv$",full.names=TRUE)
  if(is.null(SPSEL)&&!file.exists(F7r)&&length(PARTS)){write.csv(do.call(rbind,lapply(PARTS,read.csv,stringsAsFactors=FALSE)),F7r,row.names=FALSE)}
  sp<-read.csv(F7s,stringsAsFactors=FALSE);mv<-if(file.exists(F7r))read.csv(F7r,stringsAsFactors=FALSE) else NULL
  for(i in sp$spec){if(!is.null(SPSEL)&&!(i%in%SPSEL))next;if(!is.null(mv)&&i%in%mv$spec)next;vv<-strsplit(sp$vars[sp$spec==i],"\\|")[[1]];m<-safe(c(vv,FORCED),cpo=FALSE);if(is.null(m))next
    r<-do.call(rbind,lapply(vv,function(v)orr(m,v)));r$spec<-i;mv<-rbind(mv,r);if(i%%10==0){write.csv(mv,F7r,row.names=FALSE);LOG("S7 %d",i)}}
  write.csv(mv,F7r,row.names=FALSE);if(!is.null(SPSEL))quit(save="no")
  main<-read.csv(F3,stringsAsFactors=FALSE);sgn<-setNames(sign(log(main$OR)),main$cov)
  sm<-do.call(rbind,lapply(POOL,function(v){x<-mv[mv$cov==v,];ref<-sgn[v];data.frame(cov=v,label_en=LAB[v],n_specs=nrow(x),cred_frac=round(mean(x$cred),3),sign_agree=round(mean(sign(log(x$OR))==ref),3),
    OR_median=round(median(x$OR),4),OR_min=round(min(x$OR),4),OR_max=round(max(x$OR),4),robust=(mean(x$cred)>=0.75&mean(sign(log(x$OR))==ref)>=0.95))}))
  write.csv(sm,F7,row.names=FALSE);LOG("S7 robust = %s",paste(sm$cov[sm$robust],collapse=","))}

## ── S8 Knorr-Held II·III·IV(성분 인식 제약, kh_constraints.R) ──
F8<-file.path(RES,"S8_KH.csv")
if(RUN("S8")&&!file.exists(F8)){source(file.path(ROOT,"01_스크립트","kh_constraints.R"));out<-NULL
  for(ty in c("II","III","IV")){z<-kh_matrices(ty,g,11);d<-DD;d$st<-(d$time-1L)*n+d$area;stopifnot(!anyDuplicated(d$st),max(d$st)<=n*11)
    assign("KHQ",INLA::inla.as.sparse(z$Q),envir=.GlobalEnv);assign("KHC",list(A=z$A,e=rep(0,z$rankdef)),envir=.GlobalEnv)
    ex<-sprintf("f(st,model='generic0',Cmatrix=KHQ,rankdef=%d,constr=FALSE,extraconstr=KHC,hyper=pcp,values=1:%d)",z$rankdef,n*11)
    m<-timed(sprintf("S8 KH %s (rankdef %d)",ty,z$rankdef),safe(MAINV,data=d,extra=ex));if(is.null(m))next;out<-rbind(out,tab(m,MAINV,list(model=paste0("KH_",ty),WAIC=round(m$waic$waic,1))))}
  write.csv(out,F8,row.names=FALSE)}

## ── S9 지역 간·지역 내 분해(이원 중심화: 시군구 평균 = 간, 편차 = 내) ──
F9<-file.path(RES,"S9_between_within.csv")
if(RUN("S9")&&!file.exists(F9)){d<-DD;ex<-character(0);for(v in HYPV){b<-ave(d[[v]],d$area);d[[paste0(v,"__B")]]<-b;d[[paste0(v,"__W")]]<-d[[v]]-b;ex<-c(ex,bt(paste0(v,"__B")),bt(paste0(v,"__W")))}
  m<-timed("S9 between/within",safe(FORCED,data=d,extra=ex));out<-NULL;for(v in HYPV)for(k in c("B","W")){r<-orr(m,paste0(v,"__",k));r$cov<-v;r$component<-ifelse(k=="B","between","within");out<-rbind(out,r)};write.csv(out,F9,row.names=FALSE);rm(m)}

## ── S10 지연 1·2년 ──
F10<-file.path(RES,"S10_lag.csv")
if(RUN("S10")&&!file.exists(F10)){out<-NULL;for(L in 1:2){d<-DD[order(DD$area,DD$year),];for(v in HYPV)d[[v]]<-ave(d[[v]],d$area,FUN=function(s)c(rep(NA,L),head(s,-L)));d<-d[d$year>=Y0+L,];d<-d[order(d$year,d$area),];d$time<-d$year-min(d$year)+1;d$yiid<-d$time;d$id<-seq_len(nrow(d))
  m<-timed(sprintf("S10 lag %d",L),safe(MAINV,data=d));out<-rbind(out,tab(m,MAINV,list(lag=L)))};write.csv(out,F10,row.names=FALSE)}

## ── S11 시도 하나씩 빼기(그래프는 유지, 해당 시도 관측만 제외) ──
F11<-file.path(RES,"S11_leave_one_sido.csv")
if(RUN("S11")&&!file.exists(F11)){out<-NULL;for(s in sort(unique(sidoOf))){d<-DD[sidoOf[DD$area]!=s,];d$id<-seq_len(nrow(d));m<-safe(MAINV,data=d,cpo=FALSE);if(is.null(m))next;out<-rbind(out,tab(m,HYPV,list(left_out=s)))};LOG("S11 done");write.csv(out,F11,row.names=FALSE)}

## ── S12 우도 민감도 ──
F12<-file.path(RES,"S12_likelihood.csv")
if(RUN("S12")&&!file.exists(F12)){out<-NULL
  for(lk in c("binom","beta","gfree")){m<-timed(sprintf("S12 %s",lk),safe(MAINV,lik=lk,cpo=FALSE));if(is.null(m))next;out<-rbind(out,tab(m,MAINV,list(likelihood=lk)))}
  d<-DD;oa<-"eloc__DT_1YL20631__고령인구비율＜br＞(A÷B×100)";MM<-arrow::read_parquet(file.path(ROOT,"02_패널","MASTER_시군구_연도_wide.parquet"),col_select=dplyr::all_of(c("year","fbd_region",oa)))
  MM<-as.data.frame(MM);MM$fbd_region<-enc2utf8(MM$fbd_region);d$old65<-MM[[oa]][match(paste(d$year,d$fbd_region),paste(MM$year,MM$fbd_region))]
  d$old65<-ave(d$old65,d$area,FUN=function(s){k<-sum(!is.na(s));if(k==0)s else if(k==1)rep(s[!is.na(s)],length(s)) else approx(seq_along(s),s,seq_along(s),rule=2)$y});d$old65[is.na(d$old65)]<-ave(d$old65,d$year,FUN=function(s)mean(s,na.rm=TRUE))[is.na(d$old65)];stopifnot(all(!is.na(d$old65)));d$old65<-(d$old65-mean(d$old65))/sd(d$old65)
  m<-timed("S12 crude rate + 65+",safe(c(MAINV,"old65"),data=d,lik="cr",cpo=FALSE));out<-rbind(out,tab(m,MAINV,list(likelihood="crude_rate_FH_plus_old65")))
  write.csv(out,F12,row.names=FALSE)}

## ── S13 queen 그래프 · S15 PC prior 대안 ──
F13<-file.path(RES,"S13_graph_prior.csv")
if(RUN("S13")&&!file.exists(F13)){pc2<-list(prec=list(prior="pc.prec",param=c(1,0.01)),phi=list(prior="pc",param=c(0.5,0.5)))
  m1<-timed("S13 queen graph",safe(MAINV,sp=SPf("gq"),cpo=FALSE));m2<-timed("S15 PC prior U=1, phi 0.5/0.5",safe(MAINV,sp=SPf("g","pc2"),cpo=FALSE))
  write.csv(rbind(tab(m1,MAINV,list(variant="queen_graph")),tab(m2,MAINV,list(variant="pc_prior_alt"))),F13,row.names=FALSE)}

## ── S14 도농 민감도(설계 메모 §4, 적합 전 고정): 주모형 + 군 지시 + log 인구밀도 → 공간 효과·가설 OR 비교 ──
F14<-file.path(RES,"S14_rural_OR.csv")
if(RUN("S14")&&!file.exists(F14)){set.seed(4);d<-DD;d$gun<-as.numeric(grepl("군$",d$fbd_region));d$ldens<-log(d$population/(d$area_m2/1e6));d$ldens<-(d$ldens-mean(d$ldens))/sd(d$ldens)
  m0<-timed("S14 main",safe(MAINV,cpo=FALSE));m1<-timed("S14 main + county + log density",safe(c(MAINV,"gun","ldens"),data=d,cpo=FALSE))
  ir<-tab(m1,c(MAINV,"gun","ldens"));ir$variant<-"main_plus_county_density";write.csv(ir,F14,row.names=FALSE)
  r0<-spatialOR(m0);r1<-spatialOR(m1);write.csv(r1,file.path(RES,"S14_rural_spatialOR.csv"),row.names=FALSE)
  res1<-tapply(d$y_logit-(m1$summary.fitted.values$mean-m1$summary.random$id$mean[d$id]),d$area,mean);mr<-mI(as.numeric(res1),alt="two.sided")
  # 도농 R²: 시군구 평균 logit SR ~ 군 + log 밀도 + 도시 인구 비율, 잔차 Moran
  ag<-aggregate(cbind(y_logit,gun,ldens,urban_population_share)~area,d,mean);ag<-ag[order(ag$area),];lmr<-lm(y_logit~gun+ldens+urban_population_share,ag);mres<-mI(as.numeric(resid(lmr)),alt="two.sided")
  write.csv(data.frame(spearman_spatialOR_main_vs_rural=round(cor(log(r0$OR),log(r1$OR),method="spearman"),3),sd_logOR_main=round(sd(log(r0$OR)),4),sd_logOR_rural=round(sd(log(r1$OR)),4),
    hot_main=sum(r0$class=="hotspot"),hot_rural=sum(r1$class=="hotspot"),cold_main=sum(r0$class=="coldspot"),cold_rural=sum(r1$class=="coldspot"),largest_hot_rural=clusters(r1$area[r1$class=="hotspot"]),
    phi_main=fstat(m0)$phi,phi_rural=fstat(m1)$phi,WAIC_main=round(m0$waic$waic,1),WAIC_rural=round(m1$waic$waic,1),resid_moran_rural=round(mr["I"],4),resid_moran_rural_p=mr["p"],
    urbanrural_R2_pooled=round(summary(lmr)$r.squared,3),moran_after_urbanrural_lm=round(mres["I"],4),moran_after_urbanrural_lm_p=mres["p"]),file.path(RES,"S14_rural_summary.csv"),row.names=FALSE);rm(m0,m1)}
## ── S16 RSR(공간 교란, 17 B7): BYM2 총효과 ⟂ 시군구 평균 공변량 ──
F16<-file.path(RES,"S16_RSR.csv")
if(RUN("S16")&&!file.exists(F16)){Xb<-apply(as.matrix(DD[,MAINV]),2,function(x)tapply(x,DD$area,mean))[as.character(1:n),];Xb<-scale(Xb,scale=FALSE)
  RSRC<<-list(A=cbind(t(Xb),matrix(0,ncol(Xb),n)),e=rep(0,ncol(Xb)))
  m<-timed("S16 RSR",safe(MAINV,sp="f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc,extraconstr=RSRC)",cpo=FALSE));write.csv(if(is.null(m))data.frame(note="RSR fit failed (default and robust retry)") else tab(m,MAINV,list(variant="RSR")),F16,row.names=FALSE)}
## ── S17 권역 기울기(17 B6): 주모형 credible 가설 변수별 5권역 분리 ──
F17<-file.path(RES,"S17_region_slopes.csv")
if(RUN("S17")&&!file.exists(F17)&&file.exists(F3)){RG<-c("서울시"="Capital","인천시"="Capital","경기도"="Capital","대전시"="Chungcheong","세종시"="Chungcheong","충청북도"="Chungcheong","충청남도"="Chungcheong","광주시"="Honam","전라북도"="Honam","전라남도"="Honam",
    "부산시"="Yeongnam","대구시"="Yeongnam","울산시"="Yeongnam","경상북도"="Yeongnam","경상남도"="Yeongnam","강원도"="Gangwon_Jeju","제주도"="Gangwon_Jeju")
  DD$rg<-unname(RG[sidoOf[DD$area]]);stopifnot(!anyNA(DD$rg));main<-read.csv(F3,stringsAsFactors=FALSE);ROBV<-main$cov[main$cred&main$role=="hypothesis"]
  m0<-safe(MAINV,cpo=FALSE);w0<-m0$waic$waic;out<-NULL
  for(v in ROBV){D2<-DD;nv<-character(0);for(r in unique(DD$rg)){cn<-paste0(v,"__",r);D2[[cn]]<-D2[[v]]*(D2$rg==r);nv<-c(nv,cn)}
    m<-safe(c(setdiff(MAINV,v),nv),data=D2,cpo=FALSE);if(is.null(m))next;r<-do.call(rbind,lapply(nv,function(x)orr(m,x)));r$covariate<-v;r$region<-sub(".*__","",r$cov);r$dWAIC_vs_common<-round(m$waic$waic-w0,1);out<-rbind(out,r)}
  write.csv(if(is.null(out))data.frame(note="no credible hypothesis variable") else out,F17,row.names=FALSE)}
## ── S18 영향 관측(17 B10): 표준화 잔차 |r|>3 칸 제외 재적합 ──
F18<-file.path(RES,"S18_influence.csv")
if(RUN("S18")&&!file.exists(F18)){m<-safe(MAINV,cpo=FALSE);sr<-(DD$y_logit-m$summary.fitted.values$mean)/sqrt(DD$v_logit);k<-abs(sr)>3
  d<-DD[!k,];d$id<-seq_len(nrow(d));m2<-timed(sprintf("S18 drop %d cells |r|>3",sum(k)),safe(MAINV,data=d,cpo=FALSE));write.csv(tab(m2,MAINV,list(variant=sprintf("drop_%d_cells_absres_gt3",sum(k)))),F18,row.names=FALSE)}
## ── S19 조사 단절 민감도(설계 메모 §1): 2020-21 제외 · 2011-2016만 · 2017-2024만 ──
F19<-file.path(RES,"S19_breaks.csv")
if(RUN("S19")&&!file.exists(F19)){out<-NULL
  for(sp in list(list(n="drop_2020_2021",k=!DD$year%in%c(2020,2021)),list(n="2014_2019",k=DD$year<=2019),list(n="2020_2024",k=DD$year>=2020))){d<-DD[sp$k,];d$time<-d$year-min(d$year)+1;d$yiid<-d$time;d$id<-seq_len(nrow(d))
    m<-timed(sprintf("S19 %s",sp$n),safe(MAINV,data=d,cpo=FALSE));if(is.null(m))next;out<-rbind(out,tab(m,MAINV,list(variant=sp$n,n_obs=nrow(d))))}
  write.csv(out,F19,row.names=FALSE)}
## ── S20 격차 추이(설계 메모 §0 "격차가 줄었나"): 연도별 관측 P90/P10·변동계수 + 시군구별 선형 기울기(BYM2) 모형 → 기준 수준과 기울기 상관(수렴) ──
F20<-file.path(RES,"S20_gap_trend.csv")
if(RUN("S20")&&!file.exists(F20)){gy<-do.call(rbind,lapply(sort(unique(DD$year)),function(y){x<-DD$SR[DD$year==y];data.frame(year=y,mean_SR=round(mean(x),2),P90=round(quantile(x,.9),2),P10=round(quantile(x,.1),2),P90_P10=round(quantile(x,.9)/quantile(x,.1),3),P90_minus_P10=round(quantile(x,.9)-quantile(x,.1),2),CV=round(sd(x)/mean(x),3),SD_logit=round(sd(DD$y_logit[DD$year==y]),4))}))
  write.csv(gy,F20,row.names=FALSE)
  d<-DD;d$area2<-d$area;d$tc<-(d$time-mean(d$time))/10
  m<-timed("S20 area-specific slopes",safe(MAINV,data=d,cpo=FALSE,extra="f(area2,tc,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)"))
  if(!is.null(m)){ua<-m$summary.random$area$mean[1:n];us<-m$summary.random$area2$mean[1:n];ps<-sapply(m$marginals.random$area2[1:n],function(mg)1-inla.pmarginal(0,mg));set.seed(5);ms<-mI(us)
    sl<-data.frame(area=1:n,fbd_region=G$region,sido=sidoOf,level_logOR=round(ua,4),slope_per10y=round(us,4),P_slope_gt0=round(ps,4));write.csv(sl,file.path(RES,"S20_area_slopes.csv"),row.names=FALSE)
    hp<-m$summary.hyperpar;write.csv(data.frame(spearman_level_slope=round(cor(ua,us,method="spearman"),3),sd_slope=round(sd(us),4),n_slope_up=sum(ps>0.95),n_slope_down=sum(ps<0.05),moran_slope=round(ms["I"],4),moran_slope_p=ms["p"],
      phi_slope=round(if("Phi for area2"%in%rownames(hp))hp["Phi for area2","mean"] else NA,3),WAIC_slopes=round(m$waic$waic,1)),file.path(RES,"S20_slope_summary.csv"),row.names=FALSE)
    write.csv(tab(m,MAINV,list(variant="area_specific_slopes")),file.path(RES,"S20_slopes_OR.csv"),row.names=FALSE)};rm(m)}
## ── S44 역방향 사후 검정(결과 열람 후, 15 역방향 절차 2단계; 보고 시 post hoc 표시) ──
F44<-file.path(RES,"S44_reverse_posthoc.csv")
if(RUN("S44")&&!file.exists(F44)&&file.exists(F3)){main<-read.csv(F3,stringsAsFactors=FALSE);EXPD<-setNames(main$expected,main$cov)
  REV<-main$cov[main$cred&main$role=="hypothesis"&((main$OR>1&EXPD[main$cov]=="-")|(main$OR<1&EXPD[main$cov]=="+"))];out<-NULL
  ACC<-c("dentists_per_100k")
  for(v in REV){m1<-safe(c(v,FORCED),cpo=FALSE);r<-orr(m1,v);r$test<-"single_plus_forced";out<-rbind(out,r)
    m2<-safe(unique(c(v,setdiff(MAINV,setdiff(ACC,v)))),cpo=FALSE);r<-orr(m2,v);r$test<-"main_minus_other_access_vars";out<-rbind(out,r)
    m3<-safe(setdiff(MAINV,"urban_population_share"),cpo=FALSE);r<-orr(m3,v);r$test<-"main_minus_urban_share";out<-rbind(out,r)
    SES<-c("basic_livelihood_rate","tertiary_educated_share","older_single_household_share");m4<-safe(unique(c(v,setdiff(MAINV,setdiff(SES,v)))),cpo=FALSE);r<-orr(m4,v);r$test<-"main_minus_other_SES_vars";out<-rbind(out,r)
  }
  if(!is.null(out)){out$label_en<-LAB[out$cov];write.csv(out,F44,row.names=FALSE)};LOG("S44 reverse vars: %s",paste(REV,collapse=","))}
LOG("ALL DONE (%s)",if(nzchar(ONLY))ONLY else "all")
