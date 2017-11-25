(* ::Package:: *)

(* ::Title:: *)
(*Annotated code and results for "Fully Automated Cancer Diagnosis through a Tandem*)
(*of Classifiers for Digitized Histopathological Slides"*)


(* ::Section:: *)
(*Code to create the test sets and classifier functions*)


(* ::Subsection:: *)
(*Start with FCT+PCA method*)


(* ::Text:: *)
(*To access the images one can go to *)
(*https://storage.rcs-rds.ro/links/888ccd47-f84e-4510-b4dc-865cf69bca16*)
(*which may also be also found via*)
(*https://sites.google.com/site/imediatreat/data-sets*)
(*Hitting the "Save" button gives a zipped directory. I unzipped in a subdirectory of my home directory, named "histology_images" (no quotes, since I am on a Linux desktop). The images can now be imported to Mathematica as below.*)


(* ::Input:: *)
(*AbsoluteTiming[GImages=FileSystemMap[Import[#,"JPG"]&,FileNameJoin[{$HomeDirectory,"histology_images"}],{2}];]*)
(*GImages=Map[Values,Values[GImages]];*)


(* ::Input:: *)
(*Map[Length,GImages]*)


(* ::Text:: *)
(*We set up some lists for the test values (the images will change with random sampling, but not these values since we take the same amount from each of the four classes of image).*)


(* ::Input:: *)
(*testValueLens=Map[Round[Length[#]*1/3]&,GImages]*)
(*testValues=Apply[Join,MapIndexed[ConstantArray[#2[[1]]-1,#1]&,testValueLens]]*)


(* ::Text:: *)
(*Below are some parameters found to give good results for this method. We also preprocess by taking Fourier cosine transforms of all images and retaining only the low frequency square of components that was found to work best (4x4, in this case).*)


(* ::Input:: *)
(*keep=4;*)
(*dn=4;*)
(*dct=2;*)
(*GImagesS=GImages;*)


(* ::Input:: *)
(*AbsoluteTiming[dctvecs=Table[Flatten[FourierDCT[ImageData[GImagesS[[j,k]]]][[1;;dn,1;;dn]]],*)
(*{j,Length[GImages]},{k,Length[GImages[[j]]]}];]*)


(* ::Input::Initialization:: *)
nearestImages[dctvecs_,vals_,keep_]:=Module[{uu,ww,vv,udotw,norms},
{uu,ww,vv}=SingularValueDecomposition[dctvecs,keep];udotw=uu.ww;norms=(Sqrt[#1.#1]&)/@udotw;udotw=udotw/norms;
{Nearest[udotw->Transpose[{udotw,vals}]],vv}]

processInput[dctvecs_,vv_]:=
Module[
{tdotv,norms},
tdotv=dctvecs.vv;
norms=Map[Sqrt[#.#]&,tdotv];
tdotv=tdotv/norms;
tdotv]


(* ::Text:: *)
(*We use weights inversely proportional to distances between test vectors and neighbor vectors. This emulates "probabilities" of correct guesses.*)


guesses[nf_,tvecs_,n_]:=Module[
{nbrs,probs,probsB,bestvals},
probs=Table[
Module[{res=nf[tvecs[[j]],n],dists},
dists=1/Map[Norm[tvecs[[j]]-#]&,res[[All,1]]];
Thread[{res[[All,2]],dists/Total[dists]}]],
{j,Length[tvecs]}];
probsB=Map[Normal[GroupBy[#,First]]&,probs]/.(val_->vlist:{{val_,_}..}):>(val->Total[vlist[[All,2]]]);
probs=({0,1,2,3}/.probsB)/.Thread[{0,1,2,3}->0];
bestvals=Map[First[Ordering[#,1,Greater]-1]&,probs,{1}];
bestvals
]

correct[guess_,actual_]/;
Length[guess]==Length[actual]:=
Count[guess-actual,0]
correct[__]:=$Failed


(* ::Input::Initialization:: *)
runRecognitionTest[imageSets_,fraction_,repeat_,n_]:=Module[
{testSets,trainingSets,results,allTrainingImages,allTestImages,trainingValues,testValues,nfunc,vv,testvecs,guessed},
results=Table[
testSets=Map[RandomSample[#,Round[Length[#]*fraction]]&,imageSets];
trainingSets=Apply[Complement,Transpose[{imageSets,testSets}],{1}];
allTrainingImages=Apply[Join,trainingSets];
allTestImages=Apply[Join,testSets];
trainingValues=Apply[Join,MapIndexed[ConstantArray[#2[[1]]-1,Length[#1]]&,trainingSets]];
testValues=Apply[Join,MapIndexed[ConstantArray[#2[[1]]-1,Length[#1]]&,testSets]];
{nfunc,vv}=
nearestImages[allTrainingImages,trainingValues,keep];
testvecs=processInput[allTestImages,vv];
guessed=guesses[nfunc,testvecs,n];
Clear[nfunc,vv,testvecs,allTrainingImages,allTestImages];
correct[guessed,testValues]
,{j,repeat}];
{results,Mean[N@results]/Length[Apply[Join,testSets]]}
]


(* ::Text:: *)
(*We check how well the FCT+PCA method works on 100 runs. Each uses a random sample size, from the four categories, of 2/3 for training and 1/3 for testing.*)


(* ::Input:: *)
(*SeedRandom[11111];*)
(*AbsoluteTiming[runRecognitionTest[dctvecs,1/3,100,4]]*)


(* ::Text:: *)
(*So overall we are getting around 90% correct.*)


(* ::Subsection:: *)
(*This is the code to create 40 trials, using random samples of 2/3 training and 1/3 testing from each of the four categories*)


(* ::Text:: *)
(*We weight the FCT+PCA neighbors by their reciprocal distances from the input, rather than using fixed weights for the four closest. It was also tried with varying powers of those reciprocal weights but this seems to work as well as any such variant well.*)


(* ::Input:: *)
(*SeedRandom[11112222];*)
(*AbsoluteTiming[*)
(*allVals=Table[*)
(*Print[k];*)
(*testSetIndices=Map[RandomSample[Range[Length[#]],Round[Length[#]*1/3]]&,GImages];*)
(*testSets=Table[GImages[[j,testSetIndices[[j]]]],{j,Length[testSetIndices]}];*)
(*trainingSets=Apply[Complement,Transpose[{GImages,testSets}],{1}];*)
(*allTrainingImages=Apply[Join,trainingSets];*)
(*allTestImages=Apply[Join,testSets];trainingValues=Apply[Join,MapIndexed[ConstantArray[#2[[1]]-1,Length[#1]]&,trainingSets]];*)
(*testValues=Apply[Join,MapIndexed[ConstantArray[#2[[1]]-1,Length[#1]]&,testSets]];*)
(*methods={"RandomForest","NearestNeighbors","SupportVectorMachine","LogisticRegression","NeuralNetwork","NaiveBayes"};*)
(*cfuncs=Map[(Print[#];Print[AbsoluteTiming[cf=Classify[allTrainingImages->trainingValues, Method->#,PerformanceGoal->"Quality"];]];cf)&,methods];*)
(*resultprobsraw=({0,1,2,3}/.Table[Map[cfuncs[[j]][#,"Probabilities"]&,allTestImages],{j,Length[cfuncs]}]);*)
(*mlmethods[k]=resultprobsraw;*)
(*testDCTs=Table[dctvecs[[j,testSetIndices[[j]]]],{j,Length[testSetIndices]}];*)
(*trainingDCTs=Apply[Complement,Transpose[{dctvecs,testDCTs}],{1}];*)
(*allTrainingDCTs=Apply[Join,trainingDCTs];*)
(*allTestDCTs=Apply[Join,testDCTs];*)
(*{nfunc,vv}=*)
(*nearestImages[allTrainingDCTs,trainingValues,keep];*)
(*testvecs=processInput[allTestDCTs,vv];*)
(*nfuncresult=Table[*)
(*Module[{res=nfunc[testvecs[[j]],4],dists},*)
(*dists=1/Map[Norm[testvecs[[j]]-#]&,res[[All,1]]];*)
(*Thread[{res[[All,2]],dists/Total[dists]}]],*)
(*{j,Length[testvecs]}];*)
(*nfuncresultB=Map[Normal[GroupBy[#,First]]&,nfuncresult]/.(val_->vlist:{{val_,_}..}):>(val->Total[vlist[[All,2]]]);*)
(*nfuncprobs=({0,1,2,3}/.nfuncresultB)/.Thread[{0,1,2,3}->0];*)
(*fctpcaraw[k]=nfuncprobs;*)
(*Join[resultprobsraw,{nfuncprobs}]*)
(*,{k,40}];*)
(*]*)


(* ::Text:: *)
(*Nearly 5 hours later we have our data from 40 trials.*)


(* ::Input:: *)
(*allVals//Dimensions*)


(* ::Input:: *)
(*(*Export["~danl/all_values.m",allVals];*)*)


(* ::Input:: *)
(*allVals=Import["~danl/all_values.m"];*)


(* ::Subsection:: *)
(*Here are results on a per-method basis*)


(* ::Text:: *)
(*For each trial and method, use the highest probability to determine the category selected. Find the number of misses (incorrect categorizations), fractions of correct guesses, and averages of these fractions for each method, averaged over the 40 trials.*)


(* ::Input:: *)
(*testSetLengths=Map[Round[Length[#]*1/3]&,GImages]*)
(*testValues=Join@@Table[ConstantArray[j-1,testSetLengths[[j]]],{j,Length[testSetLengths]}]*)


(* ::Input:: *)
(*bestvals=Map[First[Ordering[#,1,Greater]-1]&,allVals,{3}];*)
(*diffs=Map[#-testValues&,bestvals,{2}];*)
(*misses=Map[Total,Clip[Abs[diffs]],{2}]*)
(*fracs=(119-misses)/119.*)
(*methodFracs=Mean[fracs]*)


(* ::Text:: *)
(*It seems that logistic regression and SVM are best overall, and in that order. After those we see NN, kNN, and FCT+PCA follow and are quite close to one another overall. Random forests and naive Bayes are notably below these.*)


(* ::Subsection:: *)
(*Find weights for the seven methods*)


(* ::Text:: *)
(*The hope, based on various experimenting, is that some "ensemble" apporach will be more reliable than any one method. It might even be the case that one or two methods tend to not overlap much in erroneous guesses with the other methods. These can be very useful if we can come up with a reasonable set of weights for the individual methods. Specifically, we'd like to have weights for the seven methods such that the sum of weight times probability gives a score for each of the four categories; the highest score amongst the four will determine the guessed value.*)


(* ::Text:: *)
(*The idea here is to take random subsets of 20 trials, use optimization to find "best" weights, minimizing a function of wrong guesses (we could use the raw count but we'll instead penalize misses that are off by more than a neighboring category). We can then assess the quality of the result by checking how well it works on the set withheld (think of it as a validation set). We'll do this a number of times, using different random subsets on which to perform the optimization. At the end, we look closely at the most promising to see which is our "best of the best" set of weight values.*)


(* ::Input:: *)
(*obj[vars:{_?NumberQ..},sample_]:=*)
(*With[{av=allVals[[sample]],tv=testValues},*)
(*Module[{plists,fr,diffs,misses},*)
(*plists=Map[vars.#&,av];*)
(*fr=Map[First[Ordering[#,1,Greater]-1]&,plists,{2}];*)
(*diffs=Map[#-tv&,fr];*)
(*misses=Map[Total,Abs[diffs]^2];*)
(*Total[misses]*)
(*]]*)


(* ::Input:: *)
(*testOnlyCorrect[vars:{_?NumberQ..},sample_]:=*)
(*With[{av=allVals[[sample]],tv=testValues},*)
(*Module[{plists,fr,diffs,misses},*)
(*plists=Map[vars.#&,av];*)
(*fr=Map[First[Ordering[#,1,Greater]-1]&,plists,{2}];*)
(*diffs=Map[#-tv&,fr];*)
(*misses=Map[Total,Clip[Abs[diffs]]];*)
(*Total[misses]*)
(*]]*)


(* ::Input:: *)
(*vars=Array[c,Length[allVals[[1]]]];*)
(*constraints=Flatten[{1/2<=Total[vars]<=2,Thread[0<=vars<=1]}];*)
(*high=5.;*)
(*initBounds=Map[{#,0,.3}&,vars];*)


(* ::Input:: *)
(*SeedRandom[1111122222];*)
(*wgttab=Table[*)
(*sample=Sort[RandomSample[Range[40],20]];*)
(*rest=Complement[Range[40],sample];*)
(*Print[AbsoluteTiming[bestDE=NMinimize[{obj[vars,sample],constraints},initBounds,Method->{"DifferentialEvolution","SearchPoints"->10,"PostProcess"->False,RandomSeed->RandomInteger[10^100]},MaxIterations->100]]];*)
(*Print[sampleunpenalized=testOnlyCorrect[vars/.bestDE[[2]],sample]];*)
(*Print[restunpenalized=testOnlyCorrect[vars/.bestDE[[2]],rest]];*)
(*Print[restpenalized=obj[vars/.bestDE[[2]],rest]];*)
(*{First[bestDE],restpenalized,sampleunpenalized,restunpenalized,Last[bestDE]}*)
(*,{4}]*)


(* ::Input:: *)
(*SeedRandom[1111122222];*)
(*biggerwgttabB=Table[*)
(*sample=Sort[RandomSample[Range[40],20]];*)
(*rest=Complement[Range[40],sample];*)
(*Print[AbsoluteTiming[bestDE=NMinimize[{obj[vars,sample],constraints},initBounds,Method->{"DifferentialEvolution","SearchPoints"->20,"PostProcess"->False,RandomSeed->RandomInteger[10^100]},MaxIterations->400]]];*)
(*Print[sampleunpenalized=testOnlyCorrect[vars/.bestDE[[2]],sample]];*)
(*Print[restunpenalized=testOnlyCorrect[vars/.bestDE[[2]],rest]];*)
(*Print[restpenalized=obj[vars/.bestDE[[2]],rest]];*)
(*{First[bestDE],restpenalized,sampleunpenalized,restunpenalized,Last[bestDE]}*)
(*,{20}]*)


(* ::Input:: *)
(*SeedRandom[1111122222];*)
(*biggestwgttab=Table[*)
(*sample=Sort[RandomSample[Range[40],20]];*)
(*rest=Complement[Range[40],sample];*)
(*Print[AbsoluteTiming[bestDE=NMinimize[{obj[vars,sample],constraints},initBounds,Method->{"DifferentialEvolution","SearchPoints"->50,"PostProcess"->False,RandomSeed->RandomInteger[10^100]},MaxIterations->1000]]];*)
(*Print[sampleunpenalized=testOnlyCorrect[vars/.bestDE[[2]],sample]];*)
(*Print[restunpenalized=testOnlyCorrect[vars/.bestDE[[2]],rest]];*)
(*Print[restpenalized=obj[vars/.bestDE[[2]],rest]];*)
(*{First[bestDE],restpenalized,sampleunpenalized,restunpenalized,Last[bestDE]}*)
(*,{20}]*)


(* ::Output:: *)
(*{{61.,78,40,52,{c[1]->0.0556695,c[2]->0.306545,c[3]->0.154352,c[4]->0.226961,c[5]->0.00297702,c[6]->0.00495426,c[7]->0.491235}},{51.,87,42,52,{c[1]->0.117443,c[2]->0.365041,c[3]->0.201859,c[4]->0.277663,c[5]->0.00948948,c[6]->0.,c[7]->0.615065}},{66.,67,45,52,{c[1]->0.,c[2]->0.345824,c[3]->0.3859,c[4]->0.17885,c[5]->0.,c[6]->0.00534106,c[7]->0.682427}},{67.,76,49,52,{c[1]->0.278096,c[2]->0.283574,c[3]->0.337519,c[4]->0.284524,c[5]->0.000863694,c[6]->0.,c[7]->0.694636}},{71.,71,53,42,{c[1]->0.228363,c[2]->0.299978,c[3]->0.166667,c[4]->0.223466,c[5]->0.00860804,c[6]->0.0038204,c[7]->0.483629}},{67.,61,46,46,{c[1]->0.0212204,c[2]->0.440106,c[3]->0.300166,c[4]->0.309533,c[5]->0.0043733,c[6]->0.0105019,c[7]->0.773327}},{59.,82,41,53,{c[1]->0.,c[2]->0.447378,c[3]->0.240225,c[4]->0.292026,c[5]->0.0114403,c[6]->0.00289803,c[7]->0.695566}},{70.,62,46,53,{c[1]->0.373185,c[2]->0.288958,c[3]->0.376229,c[4]->0.185511,c[5]->0.00512813,c[6]->0.,c[7]->0.664115}},{80.,56,53,44,{c[1]->0.0632418,c[2]->0.259959,c[3]->0.168696,c[4]->0.146238,c[5]->0.0306147,c[6]->0.00889573,c[7]->0.410543}},{71.,60,47,48,{c[1]->0.144832,c[2]->0.23434,c[3]->0.197326,c[4]->0.173858,c[5]->0.00718872,c[6]->0.,c[7]->0.453271}},*)
(*{54.,87,45,49,{c[1]->0.0199388,c[2]->0.403041,c[3]->0.207864,c[4]->0.311758,c[5]->0.0016314,c[6]->0.00284452,c[7]->0.664191}},{72.,66,45,49,{c[1]->0.264747,c[2]->0.387239,c[3]->0.24806,c[4]->0.253819,c[5]->0.028283,c[6]->0.00548323,c[7]->0.646465}},{59.,82,44,50,{c[1]->0.0167776,c[2]->0.311442,c[3]->0.141731,c[4]->0.22372,c[5]->0.0225016,c[6]->0.000787665,c[7]->0.496906}},{58.,87,43,52,{c[1]->0.246883,c[2]->0.444198,c[3]->0.181575,c[4]->0.315768,c[5]->0.0293493,c[6]->0.,c[7]->0.675175}},{67.,67,46,55,{c[1]->0.00013299,c[2]->0.313398,c[3]->0.362113,c[4]->0.20617,c[5]->0.000744694,c[6]->0.0117305,c[7]->0.655003}},{53.,80,47,50,{c[1]->0.033584,c[2]->0.254005,c[3]->0.292846,c[4]->0.136721,c[5]->0.,c[6]->0.00384119,c[7]->0.518084}},{68.,73,47,50,{c[1]->0.175848,c[2]->0.325059,c[3]->0.211555,c[4]->0.24828,c[5]->0.0237134,c[6]->0.00328351,c[7]->0.608238}},{59.,76,41,58,{c[1]->0.,c[2]->0.398017,c[3]->0.414852,c[4]->0.165034,c[5]->0.00224454,c[6]->0.000862637,c[7]->0.716791}},{67.,64,43,49,{c[1]->0.0111955,c[2]->0.432671,c[3]->0.260049,c[4]->0.250874,c[5]->0.0243871,c[6]->0.00173039,c[7]->0.67952}},{67.,67,49,49,{c[1]->0.00266004,c[2]->0.261016,c[3]->0.318539,c[4]->0.145525,c[5]->0.000713167,c[6]->0.0001219,c[7]->0.556737}}}*)


(* ::Text:: *)
(*Now see which optimized weights seem to be overall most promising. We expect it might be from the ones that give both sample and validation failure counts near one another, and neither much over 50.*)


(* ::Input:: *)
(*wgts=biggerwgttabB[[All,-1]];*)
(*stats=Table[*)
(*SeedRandom[1111];*)
(*bigrun=Table[testOnlyCorrect[vars/.wgts[[j]],RandomSample[Range[40],20]],{100}];*)
(*{MinMax[bigrun],Mean[N[bigrun]],Median[N@bigrun]},{j,Length[wgts]}]*)
(*beststats=Position[stats,{l1_,avg_,med_}/;avg<=47]*)


(* ::Text:: *)
(*We home in on the ones that seemed to be most promising. Of those, one is likely to be better than the rest (as we will see next when we run these over 1000 randomized subsamples).*)


(* ::Input:: *)
(*Table[*)
(*SeedRandom[1111];*)
(*biggerrun=Table[testOnlyCorrect[vars/.wgts[[j]],RandomSample[Range[40],20]],{1000}];*)
(*{MinMax[biggerrun],Mean[N[biggerrun]],Median[N@biggerrun]},{j,Flatten[beststats]}]*)


(* ::Input:: *)
(*okweights=biggerwgttabB[[12,-1]]*)


(* ::Text:: *)
(*We repeat on the longer run. The results are not too much different..*)


(* ::Input:: *)
(*wgts=biggestwgttab[[All,-1]];*)
(*stats=Table[*)
(*SeedRandom[1111];*)
(*bigrun=Table[testOnlyCorrect[vars/.wgts[[j]],RandomSample[Range[40],20]],{100}];*)
(*{MinMax[bigrun],Mean[N[bigrun]],Median[N@bigrun]},{j,Length[wgts]}]*)
(*beststats=Position[stats,{l1_,avg_,med_}/;avg<=46]*)


(* ::Input:: *)
(*Table[*)
(*SeedRandom[1111];*)
(*biggerrun=Table[testOnlyCorrect[vars/.wgts[[j]],RandomSample[Range[40],20]],{2000}];*)
(*{MinMax[biggerrun],Mean[N[biggerrun]],Median[N@biggerrun]},{j,Flatten[beststats]}]*)


(* ::Text:: *)
(*The first appears to be the best.*)


(* ::Input:: *)
(*bestweights=biggestwgttab[[1,-1]]*)


(* ::Input:: *)
(*(*bestweights={c[1]\[Rule]0.055669481826357385`,c[2]\[Rule]0.30654523340189954`,c[3]\[Rule]0.15435170605012527`,c[4]\[Rule]0.2269609192655407`,c[5]\[Rule]0.002977015367477979`,c[6]\[Rule]0.004954261132206535`,c[7]\[Rule]0.4912346173072503`};*)*)


(* ::Subsection:: *)
(*Further assessment of the results using "optimal" weights*)


(* ::Text:: *)
(*First we get the images and collapse them into a single list.*)


(* ::Input:: *)
(*allPictures=Apply[Join,GImages];*)


(* ::Text:: *)
(*Next we recreate the same random subsamples of test sets that we used in the run of 40 trials. Here each test set is 1/3 the size of the full category (levels 0-4 being the four categories). Notice we require the same RNG seed in order to do this.*)


(* ::Input:: *)
(*SeedRandom[11112222];*)
(*testSetIndices=*)
(*Table[Map[RandomSample[Range[Length[#]],Round[Length[#]*1/3]]&,GImages],*)
(*{40}];*)


(* ::Text:: *)
(*We now figure out what are the indices in the full image set that correspond to the tests.*)


(* ::Input:: *)
(*offsets=Prepend[Most[Accumulate[Map[Length,GImages]]],0]*)


(* ::Text:: *)
(*For each of the 40 trials, we have to add these offsets to the four respective categories of subsets in order to get indices of test sets in the full set of images.*)


(* ::Input:: *)
(*testIndices=Map[Apply[Join,#+offsets]&,testSetIndices];*)
(*testIndices//Dimensions*)


(* ::ItemNumbered:: *)
(*Total number of wrong guesses*)


(* ::ItemNumbered:: *)
(*Number incorrect in each of the 40 trials*)


(* ::ItemNumbered:: *)
(*Weighted incorrect (e.g. off by 2 counts as weight of 2) for each trial*)


(* ::ItemNumbered:: *)
(*For each trial, tallies of how many were off by 1, how many off by 2 (format:{{howfaroff,howmany}..})*)


(* ::ItemNumbered:: *)
(*For each miss by more than 1, the signed value of the miss (positive means we guessed high, negative means low)*)


(* ::ItemNumbered:: *)
(*A tally, over all trials, of the ordinal positions of missed guesses (tells us which slides might be commonly misclassified)*)


(* ::ItemNumbered:: *)
(*A tally, over all trials, of the ordinal positions of "badly" missed guesses (those missclassified by more than one category, that is, not put into a neighboring category).*)


(* ::ItemNumbered:: *)
(*Positions within each trial where misses occur (these are values between 1 and 119, that is, the number of tests in each trial)*)


(* ::ItemNumbered:: *)
(*Positions within each trial where large misses occur (that is, we are off by more than one category)*)


(* ::ItemNumbered:: *)
(*Ordinal positions of missed guesses in each of the 40 trials (these are the slide numbers (in range 1-357)*)


(* ::ItemNumbered:: *)
(*Ordinal positions of badly missed guesses in each of the 40 trials*)


(* ::Input:: *)
(*testAll[vars:{_?NumberQ..}]:=With[{av=allVals,tv=testValues},Module[{plists,fr,diffs,misses,rawmisses,rawmisscounts,missedPosns,bigmissPositions,bigMisses,commonMisses,commonBigMisses},*)
(*plists=Map[Total[vars*#]&,av];*)
(*keeplists=plists;*)
(*fr=Map[First[Ordering[#,1,Greater]-1]&,plists,{2}];*)
(*diffs=Map[#-tv&,fr];*)
(*misses=Map[Total,Clip[Abs[diffs]]];*)
(*rawmisses=Map[Total,Abs[diffs]];*)
(*rawmisscounts=Map[Tally,Abs[diffs]]/.{0,_}:>Nothing;*)
(*missedPosns=Map[Flatten[Position[#,aa_/;aa=!=0,Heads->False]]&,diffs];*)
(*bigmissPositions=Map[Flatten[Position[#,aa_/;Abs[aa]>=2,Heads->False],1]&,diffs];*)
(*bigMisses=MapThread[Part,{diffs,bigmissPositions}];*)
(*offBy3=Map[Flatten[Position[#,aa_/;Abs[aa]>=3,Heads->False],1]&,diffs];*)
(*commonMisses=MapThread[Extract[#1,Map[List[#]&,#2]]&,{testIndices,missedPosns}];*)
(*commonBigMisses=MapThread[Extract[#1,Map[List[#]&,#2]]&,{testIndices,bigmissPositions}];*)
(*hugebigmisses=MapThread[Extract[#1,Map[List[#]&,#2]]&,{testIndices,offBy3}];*)
(*{Total[misses],misses,rawmisses,rawmisscounts,bigMisses,SortBy[Tally[Flatten[commonMisses]],#[[2]]&],SortBy[Tally[Flatten[commonBigMisses]],#[[2]]&],missedPosns,bigmissPositions,commonMisses,commonBigMisses}*)
(*]]*)


(* ::Input:: *)
(*testAll[vars/.okweights]*)


(* ::Input:: *)
(*wbigrunStats=testAll[vars/.bestweights]*)


(* ::Text:: *)
(*We have one error in the off-by-three category.*)


(* ::Input:: *)
(*Flatten[offBy3]*)


(* ::Text:: *)
(*We'll take a look at the slide that failed so badly.*)


(* ::Input:: *)
(*Flatten[hugebigmisses]*)


(* ::Input:: *)
(*Rasterize[allPictures[[267]],ImageSize->200]*)


(* ::Text:: *)
(*Also note there were not many missed guesses that were off by 2.*)


(* ::Input:: *)
(*Total[Cases[wbigrunStats[[4]],{2,n_},{2}][[All,2]]]*)


(* ::Text:: *)
(*We next look at how many slides were ever guessed incorrectly over all 40 trials.*)


(* ::Input:: *)
(*wbigrunStats[[6]]//Length*)


(* ::Text:: *)
(*There are 29 slides that are guessed incorrectly at least once. Here I cull out the positions in the test slides of all that were guessed incorrectly more than three times.*)


(* ::Input:: *)
(*repeats1=Select[wbigrunStats[[6]],#[[2]]>3&]*)


(* ::Input:: *)
(*Total[repeats1[[All,2]]]*)


(* ::Text:: *)
(*Of 92 missed guesses, 55 came from just 7 slides. We'll show them here, in order. The first is from category 0, the second and third from category 1, and the others are from category 2. Possibly these might give some insight into the weaknesses of this methodology.*)


(* ::Input:: *)
(*badPix=Map[Rasterize[#,ImageSize->200]&,allPictures[[Sort[repeats1[[All,1]]]]]];*)


(* ::Input:: *)
(*GraphicsGrid[Partition[badPix,2,2,{1,1},{}]]*)


(* ::Subsection:: *)
(*Looking at probability vectors to assess possibly bad guesses*)


(* ::Text:: *)
(*We'll use the "good" weights to compute sets of probability "scores".*)


(* ::Input:: *)
(*Dimensions[probVectors=Map[(vars/.bestweights).#&,allVals]]*)


(* ::Text:: *)
(*Here is the first five (of 119) scores from the first of the 40 trials.*)


(* ::Input:: *)
(*probVectors[[1,1;;5]]*)


(* ::Text:: *)
(*Now we use the criterion that the largest score must be at least 0.75 of the total in order to trust a result.*)


(* ::Input:: *)
(*trustworthy=Map[Max[#]>.7*Total[#]&,probVectors,{2}];*)
(*Dimensions[trustworthy]*)
(*Tally[Flatten[trustworthy]]*)


(* ::Text:: *)
(*Now we have a look at which test results failed, to see whether most are in the "untrusted" category.*)


(* ::Input:: *)
(*missedPositions=wbigrunStats[[8]];untrustedPositions=Map[#[[All,2]]&,SplitBy[Position[trustworthy,False],First]];*)
(*checks=Table[Map[MatchQ[#,Apply[Alternatives,untrustedPositions[[j]]]]&,missedPositions[[j]]],{j,Length[missedPositions]}]*)
(*Tally[Flatten[checks]]*)


(* ::Text:: *)
(*So our criterion for being untrusted has 85 of the failures classified that way. Unfortunately that also means 7 results were wrong but marked as "trustworthy". We now check the extent to which this happens on the "big" misses where we are off by 2 or 3 categories.*)


(* ::Input:: *)
(*bigmissedPositions=wbigrunStats[[9]];*)
(*bigchecks=Table[Map[MatchQ[#,Apply[Alternatives,untrustedPositions[[j]]]]&,bigmissedPositions[[j]]],{j,Length[bigmissedPositions]}]*)
(*Tally[Flatten[bigchecks]]*)


(* ::Text:: *)
(*So none of the "big" misses were seen as trustworthy.*)


(* ::Text:: *)
(*Next we check which slides contain the  trustworthy misses.*)


(* ::Input:: *)
(*commonMisses=wbigrunStats[[-2]];*)
(*hardcases=MapThread[Pick,{commonMisses,Map[Not,checks,{2}]}]*)


(* ::Input:: *)
(*talliedhardcases=Sort[Tally[Flatten[hardcases]]]*)


(* ::Text:: *)
(*These are amongst those slides containing the most common misses (not too surprising I guess).*)


(* ::Input:: *)
(*trustedWrongPix=Map[Rasterize[#,ImageSize->200]&,allPictures[[talliedhardcases[[All,1]]]]];*)


(* ::Input:: *)
(*GraphicsGrid[Partition[trustedWrongPix,2,2,{1,1},{}]]*)


(* ::Subsection:: *)
(*Varying the probability thresholds when assessing bad guesses*)


(* ::Text:: *)
(*Here we vary the probability threshold from .5 to 1, in gradations of .05. We show for each value the tally of trusted vs untrusted tests (trusted=True), and the tally of incorrect and trusted vs incorrect and not trusted (trusted and incorrect = False). That is to say, each threshold gives a result of the form*)
(*{{{False,number not trusted},{True,number trusted}}, {{False,number incorrect but trusted},{True, number incorrect but not trusted}}}*)


(* ::Input:: *)
(*trustVsFail=Table[*)
(*trustworthy=Map[Max[#]>thresh*Total[#]&,probVectors,{2}]/.{True->"trusted results",False->"untrusted results"};*)
(*splits=SplitBy[Position[trustworthy,"untrusted results"],First];*)
(*untrustedPositions=Map[#[[All,2]]&,splits];*)
(*missing=Complement[Range[40],Flatten[Map[#[[All,1]]&,splits]]];*)
(*Do[untrustedPositions=Insert[untrustedPositions,{},j],{j,missing}];*)
(*checks=Table[Map[MatchQ[#,Apply[Alternatives,untrustedPositions[[j]]]]&,missedPositions[[j]]],{j,Length[untrustedPositions]}]/.{True->"untrusted failures",False->"trusted failures"};*)
(*trusted=FirstCase[ttally=Sort[Tally[Flatten[trustworthy]]],{"trusted results",n_}:>100*n,0];*)
(*badfailures=FirstCase[ftally=Sort[Tally[Flatten[checks]]],{"trusted failures",n_}:>n,0];*)
(*{thresh,trusted/4760.,badfailures,ttally,ftally}*)
(*,{thresh,.5,1,.05}]*)


(* ::Text:: *)
(*Note that at .8 threshold we have 73% trusted and only two trusted failures. 14 trusted failures show up when the threshold is .65, and there we trust 86.7% of the results. It seems that the threshold of .7, used in the previous section, is in the vicinity of a "sweet spot" (or "knee" of a Pareto front) in terms of having a high percentage overall of trusted results (which we want), with but few trusted failures (which of course are bad). An argument could also be made for preferring other thresholdsin the range 0.65-0.8, where the trusted failure counts are relatively low relative to trusted percentages.*)


(* ::Input:: *)
(*paretoFront=trustVsFail[[All,2;;3]]*)


(* ::Input:: *)
(*ListPlot[paretoFront,AxesLabel->{"percentage trusted","number incorrect but trusted"},AxesOrigin->{-10,-1},PlotStyle->{PointSize[Medium],Blue}]*)


(* ::Subsection:: *)
(*Confusion matrix*)


(* ::Input:: *)
(*ensembleConfusionData=Transpose[Map[Thread[{testValues,#}]&,Map[First[Ordering[#,1,Greater]-1]&,Map[(vars/.bestweights).#&,allVals],{-2}],{1}],{2,1}];*)
(*Dimensions[ensembleConfusionData]*)
(*mergedEnsmbleConfusionData=Apply[Join,ensembleConfusionData];*)
(*Dimensions[mergedEnsmbleConfusionData]*)
(*ensembleTallies=Apply[Rule[#1+1,#2]&,Tally[mergedEnsmbleConfusionData],{1}];*)
(*confusionmatrix=Normal[SparseArray[ensembleTallies]]*)


(* ::Input:: *)
(*tg=Labeled[TextGrid[Prepend[Join[Transpose@{Range[0,3]},confusionmatrix,2],Prepend[Range[0,3],""]],Frame->All,Dividers->{{{Thickness[2],Thickness[2],True,True,True}},{{Thickness[2],Thickness[2],True,True,True}}},Background->{1->White,{1->White},Join[{{{2,-1},{2,-1}}->LightRed},Thread[Map[{#,#}&,Range[2,5]]->Green]]}],{"expected","actual","confusion matrix"},{Left,Top,Bottom}]*)


(* ::Subsection:: *)
(*Apply the optimal weights to each of the 40 individual trials*)


(* ::Text:: *)
(*We show the number of misses, and correctness percentage, for each of the 40 trials, as assessed using the best set of weights. We then give the range, mean, and median of this set of correct percentages.*)


(* ::Input:: *)
(*totaledValues=Map[(vars/.bestweights).#&,allVals];*)
(*bestvals=Map[First[Ordering[#,1,Greater]-1]&,totaledValues,{2}];*)
(*diffs=Map[#-testValues&,bestvals,{1}];*)
(*misses=Map[Total,Clip[Abs[diffs]],{1}];*)
(*correctPercentages=(119-misses)/119.;*)
(*Transpose[{misses,correctPercentages}]*)
(*{MinMax[correctPercentages],Mean[correctPercentages],Median[correctPercentages]}*)


(* ::Text:: *)
(*Perhaps notable is that even the worst trial still came to 94% correct.*)


(* ::Subsection:: *)
(*Other things attempted*)


(* ::Text:: *)
(*I tried to forego optimizing weights, instead using reciprocals of individual method success rates (as measured over the 40 trials). Variant: use some power of those reciprocals (helps to separate them better).*)
(*Also used this in tandem with raising probabilities to some power, such as 2. The idea here is to better emphasize results that seem "certain" over those that are guessing lower probabilities.*)
(*Upshot: I found no combination of these that gave better than a 97% success rate overall. So the DE-based weight computation seems to be important.*)


(* ::Subsection:: *)
(*Correlations between methods*)


(* ::Input:: *)
(*allVals=Import["all_values.m"];*)
(*tvals=Transpose[allVals,{2,1}];*)
(*Dimensions[tvals]*)
(*tjvals=Apply[Join,tvals,{1}];*)
(*Dimensions[tjvals]*)


(* ::Input:: *)
(*correlations=Table[{i,j,Correlation[tjvals[[i]],tjvals[[j]]]},{i,Length[tjvals]},{j,i}]*)


(* ::Text:: *)
(*Now we flatten and check correlations.*)


(* ::Input:: *)
(*tjfvals=Map[Flatten,tvals];*)
(*Dimensions[tjfvals]*)


(* ::Input:: *)
(*correlations=Table[{i,j,Correlation[tjfvals[[i]],tjfvals[[j]]]},{i,Length[tjvals]},{j,i}]*)


(* ::Subsubsection:: *)
(*Cohen Kappa correlation (unweighted)*)


(* ::Text:: *)
(*Here we show the Cohen Kappa correlation values for all pairs of methods.*)


(* ::Input:: *)
(*totaledValues=Map[(vars/.naiveweights).#&,allVals];*)
(*bestvals7=Map[First[Ordering[#,1,Greater]-1]&,allVals,{3}];*)


(* ::Input:: *)
(*bestvalsLists=Transpose[Flatten[Transpose[bestvals7,{1,3,2}],1]];*)
(*Dimensions[bestvalsLists]*)


(* ::Input:: *)
(*cohenKappa[l1_,l2_,rng_]:=Module[*)
(*{tallies=Tally[Transpose[{l1,l2}]],tally1Vec,tally2Vec,len=Length[l1],agreed,pObs,pEmp},*)
(*agreed=Total[Cases[tallies,{{a_,a_},b_}:>b]];*)
(*pObs=agreed/len;*)
(*tally1Vec=SparseArray[Map[#[[1]]+1-rng[[1]]->#[[2]]&,Tally[l1]],Length[rng]];*)
(*tally2Vec=SparseArray[Map[#[[1]]+1-rng[[1]]->#[[2]]&,Tally[l2]],Length[rng]];*)
(*pEmp=tally1Vec.tally2Vec/len^2;*)
(*(pObs-pEmp)/(1-pEmp)*)
(*]*)


(* ::Input:: *)
(*cohenKappaTable=Table[{i,j,cohenKappa[bestvalsLists[[i]],bestvalsLists[[j]],Range[0,3]]},{i,Length[bestvalsLists]},{j,i}];*)
(*N[cohenKappaTable]*)
