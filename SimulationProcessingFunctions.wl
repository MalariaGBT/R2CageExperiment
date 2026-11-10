(* ::Package:: *)

SetDirectory["~/Dropbox/Doublesex/qfs5/SouthKen/Nov26"];
(*-----------Import experimental data--------------------*)
<<"Pupadata";
<<"ExPlots";
(*--------------------------------------------------------------------------------------------*)
(*Function that takes simulation output 'res' and samples from it phenotypic data corresponding to the empirical data.
 sampsize is the number of eggs set aside for screening when they are developed into pupae, and larvmort is larval mortality*)
phenotype[res_,sampsize_,larvmort_]:=Block[{ncage,ngen,nsamps,phen,phentake,tots,fem,rfp,gfp},
ncage=Length[res];
ngen=Length[res[[1]]];
nsamps:=Random[BinomialDistribution[sampsize,1-larvmort]];
phen=Table[Map[Total,{res[[cage,gen,{6,12}]],res[[cage,gen,{5,10,11}]],res[[cage,gen,{3,8,9}]],res[[cage,gen,{7}]],res[[cage,gen,{4}]],res[[cage,gen,{2}]],res[[cage,gen,{1}]]}],{cage,ncage},{gen,ngen}];
phentake=Table[RandomVariate[MultinomialDistribution[nsamps,phen[[cage,gen]]/Total[phen[[cage,gen]]]]],{cage,ncage},{gen,ngen}];
tots=Table[Map[Total,phentake[[cage]]],{cage,ncage}];
fem=Table[1. Total[phentake[[cage,gen,{5,6,7}]]]/tots[[cage,gen]],{cage,ncage},{gen,ngen}];
rfp=Table[1. Total[phentake[[cage,gen,{1,3,6}]]]/tots[[cage,gen]],{cage,ncage},{gen,ngen}];
gfp=Table[1. Total[phentake[[cage,gen,{1,2,5}]]]/tots[[cage,gen]],{cage,ncage},{gen,ngen}];
{phentake,fem,rfp,gfp}];
(*--------------------------------------------------------------------------------------------*)

(*--------------------------------------------------------------------------------------------*)
(*Function to compute sum of squares distance from simulated to experimental data*)
Dis[SimData_,CageData_]:=Sum[(SimData[[cage,gen,phen]]-CageData[[cage,gen,phen]])^2,{cage,Length[SimData]},{gen,2,Length[SimData[[1]]],1},{phen,7}];
(*--------------------------------------------------------------------------------------------*)
(*--------------------------data variance metric, used in post-predictive checks--------------*)
T[da_]:=Block[{timeMean},
timeMean=Table[N[Mean[da[[cage,2;;Dimensions[da][[2]],phen]]]],{cage,Dimensions[da][[1]]},{phen,Dimensions[da][[3]]}];
Sum[(da[[cage,time,phen]]-timeMean[[cage,phen]])^2,{time,2,Dimensions[da][[2]],1},{cage,Dimensions[da][[1]]},{phen,Dimensions[da][[3]]}]];
(*--------------------------------------------------------------------------------------------*)

(*---Function that samples the prior, generates simulated data, and computes the distance to the experimental data-------*)
GenerateParticle[pars_,pupae_,index_]:=Block[{parTemp=pars,cpp,LarvMort,pp,fem,rfp,gfp,Tsim},
parTemp[[2]]=.5Random[];
LarvMort=0.5Random[];
Export["./Par"<>ToString[index]<>".csv",parTemp,"CSV","TextDelimiters"->""];
cpp=Partition[Import["!./CageSimulation<Par"<>ToString[index]<>".csv","Table"],pars[[9]]+1];

{pp,fem,rfp,gfp}=phenotype[cpp,200,LarvMort];
Tsim=T[pp];
{parTemp[[2]],LarvMort,Dis[pp,pupae],Tsim}];
(*--------------------------------------------------------------------------------------------*)
(*--Function that generates plots of simulations Vs data, by taking 'NumDraws' samples from the Posterior distributin 'post'-----*)
GenerateDataPlot[pars_,post_,NumDraws_,index_,th_]:=Block[{parSamp=pars,mu,cpp,pp,fem,rfp,gfp,femAll={},rfpAll={},gfpAll={},
plots,femplotsAll,rfpplotsAll,gfpplotsAll},
	
Do[
{parSamp[[2]],mu}=RandomSample[post[[All,{1,2}]],1][[1]];
Export["./Par"<>ToString[index]<>".csv",parSamp,"CSV","TextDelimiters"->""];cpp=Partition[Import["!../cpp/testB<Par"<>ToString[index]<>".csv","Table"],pars[[9]]+1];
{pp,fem,rfp,gfp}=phenotype[cpp,200,mu];
AppendTo[femAll,fem];AppendTo[rfpAll,rfp];AppendTo[gfpAll,gfp],{ii,NumDraws}];
	plots=Table[Show[Table[ListLinePlot[{femAll[[samp,cage]],rfpAll[[samp,cage]],gfpAll[[samp,cage]]}[[phen]],
	DataRange->{0,8},PlotStyle->{{(*Lighter[Blue,.7Random[]],*)GrayLevel[0.25+0.5Random[]],Thickness[th]},{(*Lighter[Red,.7Random[]],*)
	GrayLevel[0.25+0.5Random[]],Thickness[th]},{(*Lighter[Green,.7Random[]],*)
	GrayLevel[0.25+0.5Random[]],Thickness[th]}}[[phen]]],{samp,NumDraws}]],{phen,3},{cage,Length[femAll[[1]]]}];
femplotsAll=Table[Show[{plots[[1,cage]],femplots[[cage]]},PlotRange->{{-.5,8.5},{0,1}},
AspectRatio->1.5,Frame->True,PlotStyle->Black,Axes->None,FrameLabel->{"Generation",None}],{cage,3}];
rfpplotsAll=Table[Show[{plots[[2,cage]],rfpplots[[cage]]},PlotRange->{{-.5,8.5},{0,1}},
AspectRatio->1.5,Frame->True,PlotStyle->Black,Axes->None,FrameLabel->{"Generation",None}],{cage,3}];
gfpplotsAll=Table[Show[{plots[[3,cage]],gfpplots[[cage]]},PlotRange->{{-.5,8.5},{0,1}},
AspectRatio->1.5,Frame->True,PlotStyle->Black,Axes->None,FrameLabel->{"Generation",None}],{cage,3}];
GraphicsGrid[{femplotsAll,rfpplotsAll,gfpplotsAll}]];
(*--------------------------------------------------------------------------------------------*)
