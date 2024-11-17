## INPUT(S)
# INCHIKEYs_CNL_Ref_PestMix_1-8.csv
# PestMix1-8_test_report_comp_IDs.csv
# dataAllFP_withNewPredictedRiWithStratification.csv
# TPTN_dfCNLfeaturesStr.csv
# CocamideExtended73_CNLsRi_RFwithStratification.joblib

## OUTPUT(S)
# PestMix1-8_test_report_comp_IDs_check.csv
# PestMix1-8_test_report_comp_IDs_checked.csv
# PestMix1-8_test_report_comp_IDs_ready4CNLdf.csv
# PestMix1-8_test_report_comp_IDs_extractedWithCNLsList.csv
# PestMix1-8_test_report_comp_IDs_withCNLRideltaRi.csv
# PestMix1-8_test_report_comp_IDs_dataframeTPTNModeling.csv

VERSION
## install packages needed ##
using Pkg
#Pkg.add("ScikitLearn")
#Pkg.add(PackageSpec(url=""))

## import packages from Julia ##
import Conda
Conda.PYTHONDIR
ENV["PYTHON"] = raw"C:\Users\user\AppData\Local\Programs\Python\Python311\python.exe"  # python 3.11
Pkg.build("PyCall")
Pkg.status()
using Random
using CSV, DataFrames, Conda, LinearAlgebra, Statistics
using PyCall
using StatsPlots
using Plots
using ProgressBars
using ScikitLearn  #: @sk_import, fit!, predict
@sk_import ensemble: RandomForestRegressor
@sk_import ensemble: GradientBoostingClassifier
@sk_import linear_model: LogisticRegression
@sk_import ensemble: RandomForestClassifier
@sk_import ensemble: AdaBoostClassifier
@sk_import tree: DecisionTreeClassifier
@sk_import metrics: recall_score
@sk_import neural_network: MLPClassifier
@sk_import svm: SVC
@sk_import neighbors: KNeighborsClassifier
@sk_import inspection: permutation_importance
#using ScikitLearn.GridSearch: RandomizedSearchCV
using ScikitLearn.CrossValidation: cross_val_score
using ScikitLearn.CrossValidation: train_test_split
#using ScikitLearn.GridSearch: GridSearchCV

## import packages from Python ##
jl = pyimport("joblib")             # used for loading models
f1_score = pyimport("sklearn.metrics").f1_score
matthews_corrcoef = pyimport("sklearn.metrics").matthews_corrcoef
make_scorer = pyimport("sklearn.metrics").make_scorer
f1 = make_scorer(f1_score, pos_label=1, average="binary")

## import groud truth ##
INCHIKEYreal = Array(CSV.read("F:\\UvA\\INCHIKEYs_CNL_Ref_PestMix_1-8.csv", DataFrame)[:,1])

## handle MS/MS data ##
    ## input __ x 20 dfs -> 29914 x 12+4+1 df
    # ID, Rt, MS1Mass, Name, Formula, ACCESSION, 
    # RefMatchFrag, UsrMatchFrag, MS1Error, MS2Error, MS2ErrorStd, 
    # DirectMatch, ReversMatch, Probability, FinalScore, 
    # SpecType, MatchedFrags, Inchikey, FragMZ, FragInt
    inputDB1 = CSV.read("F:\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs.csv", DataFrame)
    #inputDB5 = CSV.read("F:\\Cand_synth_rr10_4001_5000.csv", DataFrame)
    #combinedDB = vcat(inputDB1, inputDB2, inputDB3, inputDB4, inputDB5)
    combinedDB = inputDB1[:, ["ID", 
        #"INCHIKEYreal", "INCHIKEY", 
        "Inchikey", 
        "RefMatchFrag", "UsrMatchFrag", 
        "MS1Error", "MS2Error", "MS2ErrorStd", 
        "DirectMatch", "ReversMatch", 
        "FinalScore", "MS1Mass", "FragMZ"]]
    combinedDB = combinedDB[combinedDB.RefMatchFrag .!= "NaN", :]
    combinedDB = combinedDB[combinedDB.RefMatchFrag .!== NaN, :]
    combinedDB = combinedDB[combinedDB.UsrMatchFrag .!= "NaN", :]
    combinedDB = combinedDB[combinedDB.UsrMatchFrag .!== NaN, :]
    combinedDB = combinedDB[combinedDB.MS1Error .!= "NaN", :]
    combinedDB = combinedDB[combinedDB.MS1Error .!== NaN, :]
    combinedDB = combinedDB[combinedDB.MS2Error .!= "NaN", :]
    combinedDB = combinedDB[combinedDB.MS2ErrorStd .!== NaN, :]
    combinedDB = combinedDB[combinedDB.DirectMatch .!= "NaN", :]
    combinedDB = combinedDB[combinedDB.DirectMatch .!== NaN, :]
    combinedDB = combinedDB[combinedDB.ReversMatch .!= "NaN", :]
    combinedDB = combinedDB[combinedDB.ReversMatch .!== NaN, :]
    combinedDB = combinedDB[combinedDB.FinalScore .!= "NaN", :]
    combinedDB = combinedDB[combinedDB.FinalScore .!== NaN, :]
    combinedDB = combinedDB[combinedDB.Inchikey .!= "N/A", :]
    combinedDB[!, "RefMatchFragRatio"] .= float(0)
    combinedDB[!, "UsrMatchFragRatio"] .= float(0)
    combinedDB[!, "FinalScoreRatio"] .= float(0)
    combinedDB[!, "INCHIKEY_ID"] .= ""
    
## handle 8 features ##
    function takeRatio(str)
        num = ""
        ratio = float(0)
        for i in 1:length(str)
            if num == "-"
                num = ""
                continue
            elseif str[i] == '-'
                ratio = parse(Float64, num)
                #print(ratio)
                num = "-"
            elseif str[i] != '-'
                num = string(num, str[i])
                #print(num)
            end
        end
        ratio = ratio / (parse(Float64, num))
        return ratio
    end
    #
    ratioRef = []
    ratioUsr = []
    ratioScore = []
    trueOrFalse = []
    for i in 1:size(combinedDB, 1)
        println(i)
        push!(ratioRef, takeRatio(combinedDB[i, "RefMatchFrag"]))
        push!(ratioUsr, takeRatio(combinedDB[i, "UsrMatchFrag"]))
        push!(ratioScore, combinedDB[i, "FinalScore"]/7)
        combinedDB[i, "INCHIKEY_ID"] = string(combinedDB[i, "Inchikey"], "_", string(combinedDB[i, "ID"]))
        if (String31(combinedDB[i, "Inchikey"]) in INCHIKEYreal)
            push!(trueOrFalse, 1)
        else
            push!(trueOrFalse, 0)
        end
    end
    #
    combinedDB[!, "RefMatchFragRatio"] = ratioRef
    combinedDB[!, "UsrMatchFragRatio"] = ratioUsr
    combinedDB[!, "FinalScoreRatio"] = ratioScore
    combinedDB[!, "LABEL"] = trueOrFalse
    savePath = "F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_check.csv"
    CSV.write(savePath, combinedDB)
    #
    outputDf = combinedDB[:, ["INCHIKEY_ID", "Inchikey", "ID", "RefMatchFragRatio", "UsrMatchFragRatio", 
    "MS1Error", "MS2Error", "MS2ErrorStd", "DirectMatch", "ReversMatch", 
    "FinalScoreRatio", "MS1Mass", "FragMZ", "LABEL"]]
## output csv ## is a __ x 14 df
savePath = "F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_checked.csv"
CSV.write(savePath, outputDf)

## handle the delta Ri feature ##
    ## create a __ x 2+8+2+1+1 df ##
        # columns: INCHIKEY_ID, INCHIKEYreal, 8+1 ULSA features, LABEL
        #                      FP->Ri, CNL->Ri ^
    ## find FP-derived Ri values ##
    ## match INCHIKEY ##, 30684 x 793 df
    inputFP2Ri = CSV.read("F:\\UvA\\dataAllFP_withNewPredictedRiWithStratification.csv", DataFrame)
    sort!(inputFP2Ri, [:INCHIKEY, :SMILES])
    #
    outputDf[!, "FPpredictRi"] .= float(0)
    for i in 1:size(outputDf, 1)
        println(i)
        if outputDf[i, "Inchikey"] in Array(inputFP2Ri[:, "INCHIKEY"])
            rowNo = findall(inputFP2Ri.INCHIKEY .== outputDf[i, "Inchikey"])
            outputDf[i, "FPpredictRi"] = inputFP2Ri[rowNo[end:end], "predictRi"][1]
        else
            outputDf[i, "FPpredictRi"] = float(8888888)
        end
    end
    ## filter in INCHIKEY_ID with Ri values ##
    outputDf = outputDf[outputDf.FPpredictRi .!= float(8888888), :]
    sort!(outputDf, [:LABEL, :INCHIKEY_ID])
## output csv __ x 15 df ##
savePath = "F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_ready4CNLdf.csv"
CSV.write(savePath, outputDf)

## create CNL df ##
inputTPTNdf = CSV.read("F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_ready4CNLdf.csv", DataFrame)
    ## define a function for data extraction ##
    function getVec(matStr)
        if matStr[1] .== '['
            if contains(matStr, ", ")
                str = split(matStr[2:end-1],", ")
            else
                str = split(matStr[2:end-1]," ")
            end
        elseif matStr[1] .== 'A'
            if contains(matStr, ", ")
                str = split(matStr[5:end-1],", ")
            else
                str = split(matStr[5:end-1]," ")
            end
        elseif matStr[1] .== 'F'
            if matStr .== "Float64[]"
                return []
            else
                str = split(matStr[9:end-1],", ")
            end
        elseif matStr[1] .== 'I'
            if matStr .== "Int64[]"
                return []
            else
                str = split(matStr[7:end-1],", ")
            end
        else
            println("New vector start")
            println(matStr)
        end
        if length(str) .== 1 && cmp(str[1],"") .== 0
            return []
        else
            str = parse.(Float64, str)
            return str
        end
    end

    ## initialize an aray for 1 more column ## -> __ x 15+1
    inputTPTNdf[!, "CNLmasses"] .= [[]]
    size(inputTPTNdf)

    ## import pre-defined 15961 CNL features ##
    CNLfeaturesStr = CSV.read("F:\\TPTN_dfCNLfeaturesStr.csv", DataFrame)[:, "CNLfeaturesStr"]
        #
        candidatesList = []
        for can in CNLfeaturesStr
            #push!(candidatesList, round(parse(Float64, can), digits = 2))
            push!(candidatesList, round(can, digits = 2))
        end
        CNLfeaturesStr = []
        for can in candidatesList
            #push!(candidatesList, round(parse(Float64, can), digits = 2))
            push!(CNLfeaturesStr, string(can))
        end

    ## extract CNL ##
    for i in 1:size(inputTPTNdf, 1)
        println(i)
        fragIons = getVec(inputTPTNdf[i,"FragMZ"])
        arrNL = Set()
        for frag in fragIons
            if (inputTPTNdf[i,"MS1Mass"] - frag) >= float(0)
                NL = round((inputTPTNdf[i,"MS1Mass"] - frag), digits = 2)
                if (NL in candidatesList)
                    push!(arrNL, NL)
                end
            end
        end
        inputTPTNdf[i, "CNLmasses"] = sort!(collect(arrNL))
    end
    sort!(inputTPTNdf, [:LABEL, :INCHIKEY_ID, :CNLmasses])
    
    ## reduce df size (rows) ##
    function getMasses(db, i, arr, arrType = "dig")
        massesArr = arr
        if (arrType == "dig")
            masses = db[i, "CNLmasses"]
        elseif (arrType == "str")
            masses = getVec(db[i, "CNLmasses"])
        end
        for mass in masses
            push!(massesArr, mass)
        end
        return massesArr
    end

    ## remove rows that has Frag-ion of interest < 2 (optional) ##
    retain = []
    for i in 1:size(inputTPTNdf, 1)
        println(i)
        arr = []
        arr = getMasses(inputTPTNdf, i, arr, "dig")#"str")
        if (size(arr, 1) >= 2)
            push!(retain, i)
        end
    end
    inputTPTNdf = inputTPTNdf[retain, :]
    savePath = "F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_extractedWithCNLsList.csv"
    CSV.write(savePath, inputTPTNdf)

    ## store data in a Matrix ##
    inputTPTNdf = CSV.read("F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_extractedWithCNLsList.csv", DataFrame)
    X = zeros(18183, 15961)
    #
    for i in 1:18183
        println(i)
        arr = []
        arr = getMasses(inputTPTNdf, i, arr, "str")#"dig")
        mumIon = round(inputTPTNdf[i, "MS1Mass"], digits = 2)
        for col in arr
            mz = findall(x->x==col, candidatesList)
            if (col <= mumIon)
                X[i, mz] .= 1
            elseif (col > mumIon)
                X[i, mz] .= -1
            end
        end
    end

    ## create df ## with 11 + 1(monoisotopic mass) + 15961(CNLs)
    dfCNLs = DataFrame(X, CNLfeaturesStr)
        insertcols!(dfCNLs, 1, ("ENTRY"=>collect(1:18183)))
        insertcols!(dfCNLs, 2, ("INCHIKEY_ID"=>inputTPTNdf[:, "INCHIKEY_ID"]))
        insertcols!(dfCNLs, 3, ("INCHIKEY"=>inputTPTNdf[:, "Inchikey"]))
        insertcols!(dfCNLs, 4, ("INCHIKEYreal"=>inputTPTNdf[:, "ID"]))
        insertcols!(dfCNLs, 5, ("RefMatchFragRatio"=>inputTPTNdf[:, "RefMatchFragRatio"]))
        insertcols!(dfCNLs, 6, ("UsrMatchFragRatio"=>inputTPTNdf[:, "UsrMatchFragRatio"]))
        insertcols!(dfCNLs, 7, ("MS1Error"=>inputTPTNdf[:, "MS1Error"]))
        insertcols!(dfCNLs, 8, ("MS2Error"=>inputTPTNdf[:, "MS2Error"]))
        insertcols!(dfCNLs, 9, ("MS2ErrorStd"=>inputTPTNdf[:, "MS2ErrorStd"]))
        insertcols!(dfCNLs, 10, ("DirectMatch"=>inputTPTNdf[:, "DirectMatch"]))
        insertcols!(dfCNLs, 11, ("ReversMatch"=>inputTPTNdf[:, "ReversMatch"]))
        insertcols!(dfCNLs, 12, ("FinalScoreRatio"=>inputTPTNdf[:, "FinalScoreRatio"]))
        insertcols!(dfCNLs, 13, ("MONOISOTOPICMASS"=>((inputTPTNdf[:, "MS1Mass"] .- 1.007276)/1000)))
    dfCNLs[!, "FPpredictRi"] = inputTPTNdf[:,"FPpredictRi"]
    size(dfCNLs)  # __ x (13+15961+1)

    ## Load a pre-trained CNL-to-Ri model ##
    # requires python 3.11 or 3.12
    modelRF_CNL = jl.load("F:\\UvA\\CocamideExtended73_CNLsRi_RFwithStratification.joblib")
    size(modelRF_CNL)
    
    ## predict CNL-derived Ri ##
    CNLpredictedRi = predict(modelRF_CNL, Matrix(dfCNLs[:, 13:end-1]))
    dfCNLs[!, "CNLpredictRi"] = CNLpredictedRi
    dfCNLs[!, "DeltaRi"] = (CNLpredictedRi - dfCNLs[:, "FPpredictRi"]) / 1000
    dfCNLs[!, "LABEL"] = inputTPTNdf[:, "LABEL"]
## save csv ##
savePath = "F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_withCNLRideltaRi.csv"
CSV.write(savePath, dfCNLs)
println("done for saving csv")

## apply leverage filter ##
dfOutput = CSV.read("F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_withCNLRideltaRi.csv", DataFrame)
    #
    describe(dfOutput)[end-4:end-2, :]
    X = deepcopy(dfOutput[:, 13:end-4])  # __ x 15962 df
    size(X)
    Y = deepcopy(dfOutput[:, end-2])
    size(Y)
    Xmat = Matrix(X)
    #
    # 15962 x 15962
    hipinv = zeros(15962, 15962)
    hipinv[:,:] .= pinv(Xmat'*Xmat)
    #
    function leverage_dist(X)
        h = zeros(18183,1)
        for i in ProgressBar(1: size(X,1)) #check dimensions
            x = X[i,:]
            #hi = x'*pinv(X'*X)*x
            hi = x'* hipinv *x
            #push!(h,hi)
            h[i,1] = hi
        end
        return h
    end
    #
    h = leverage_dist(Matrix(X))
    ht = Vector(transpose(h)[1,:])
    dfOutput[!, "Leverage"] .= ht
    dfOutput = dfOutput[:, vcat(collect(1:13), end-4, end-3, end-2, end-1, end)]
## save csv ##
savePath = "F:\\UvA\\app\\PestMix1-8_1000ug-L_NoTea_1ul_AllIon_pos_52_report_comp_IDs_dataframeTPTNModeling.csv"
CSV.write(savePath, dfOutput)
