## INPUT(S)
# PestMix1-8_test_report_comp_IDs_dataframeTPTNModeling.csv
# modelTPTNModeling_6paraKNN_noFilterWithDeltaRI.joblib

## OUTPUT(S)
# PestMix1-8_test_report_comp_IDs_dataframeTPTNModeling_0d5FinalScoreRatioDEFilterSTD.csv
# PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTN_KNN.csv
# PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTNandpTP_KNN.csv
# PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTNandpTP_KNN_ind.csv

VERSION
## install packages needed ##
using Pkg
#Pkg.add("ScikitLearn")
#Pkg.add(PackageSpec(url=""))

## import packages from Julia ##
import Conda
Conda.PYTHONDIR
ENV["PYTHON"] = raw"C:\Users\T1208\AppData\Local\Programs\Python\Python311\python.exe"  # python 3.11
Pkg.build("PyCall")
Pkg.status()
using Random
using CSV, DataFrames, Conda, LinearAlgebra, Statistics
using PyCall
using StatsPlots
using Plots
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

## input and pre-process sample TP/TN set ##
inputDB_test = CSV.read("F:\\UvA\\app\\PestMix1-8_test_report_comp_IDs_dataframeTPTNModeling.csv", DataFrame)
    sort!(inputDB_test, [:ENTRY])
    insertcols!(inputDB_test, 10, ("MatchDiff"=>float(0)))
    inputDB_test = inputDB_test[inputDB_test.FinalScoreRatio .>= float(0.5), :]
    inputDB_test = inputDB_test[inputDB_test.Leverage .<= 0.14604417882015916, :]
    for i = 1:size(inputDB_test, 1)
        inputDB_test[i, "RefMatchFragRatio"] = log10(inputDB_test[i, "RefMatchFragRatio"])
        inputDB_test[i, "UsrMatchFragRatio"] = log10(inputDB_test[i, "UsrMatchFragRatio"])
        inputDB_test[i, "FinalScoreRatio"] = log10(inputDB_test[i, "FinalScoreRatio"])
        inputDB_test[i, "MatchDiff"] = inputDB_test[i, "DirectMatch"] - inputDB_test[i, "ReversMatch"]
        inputDB_test[i, "MONOISOTOPICMASS"] = log10(inputDB_test[i, "MONOISOTOPICMASS"])
        if inputDB_test[i, "DeltaRi"] !== float(0)
            inputDB_test[i, "DeltaRi"] = inputDB_test[i, "DeltaRi"] * float(-1)
        end
    end
    #
    for f = 5:14
        avg = float(mean(inputDB_test[:, f]))
        top = float(maximum(inputDB_test[:, f]))
        down = float(minimum(inputDB_test[:, f]))
        for i = 1:size(inputDB_test, 1)
            inputDB_test[i, f] = (inputDB_test[i, f] - avg) / (top - down)
        end
    end

## save  ##, ouputing __ x 18 df
savePath = "F:\\UvA\\app\\PestMix1-8_test_report_comp_IDs_dataframeTPTNModeling_0d5FinalScoreRatioDEFilterSTD.csv"
CSV.write(savePath, inputDB_test)
inputDB_test[inputDB_test.LABEL .== 1, :]


# ==================================================================================================
## define functions for performace evaluation ##
    # Maximum absolute error
    # mean square error (MSE) calculation
    # Root mean square error (RMSE) calculation
    function errorDetermination(arrRi, predictedRi)
        sumAE = 0
        maxAE = 0
        for i = 1:size(predictedRi, 1)
            AE = abs(arrRi[i] - predictedRi[i])
            if (AE > maxAE)
                maxAE = AE
            end
            sumAE += (AE ^ 2)
        end
        MSE = sumAE / size(predictedRi, 1)
        RMSE = MSE ^ 0.5
        return maxAE, MSE, RMSE
    end
    #
    # Average score
    function avgScore(arrAcc, cv)
        sumAcc = 0
        for acc in arrAcc
            sumAcc += acc
        end
        return sumAcc / cv
    end


# ==================================================================================================
## deploy models for training set ##
    ## load a model ##, with DeltaRi
    # requires python 3.11 or 3.12
    model = jl.load("F:\\UvA\\modelTPTNModeling_6paraKNN_noFilterWithDeltaRI.joblib")
    #
    ## deploy model ##
    predictedTPTN_test = predict(model, Matrix(inputDB_test[:, vcat(5, 7,9, 13,14, 17)]))
    inputDB_test[!, "withDeltaRIpredictTPTN"] = predictedTPTN_test
    #
## save ##, ouputing testSet df __ x 19 df
savePath = "F:\\UvA\\app\\PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTN_KNN.csv"
CSV.write(savePath, inputDB_test)

## evaluate predictive performance of the model with DeltaRi ##
inputTestDB_withDeltaRiTPTN = CSV.read("F:\\UvA\\app\\PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTN_KNN.csv", DataFrame)
    #
    maxAE_val, MSE_val, RMSE_val = errorDetermination(inputTestDB_withDeltaRiTPTN[:, end-2], inputTestDB_withDeltaRiTPTN[:, end])
    pTP_test = predict_proba(model, Matrix(inputTestDB_withDeltaRiTPTN[:, vcat(5, 7,9, 13,14, 17)]))  # __ × 2 Matrix
    recall_test = recall_score(Vector(inputTestDB_withDeltaRiTPTN[:, end-2]), predict(model, Matrix(inputTestDB_withDeltaRiTPTN[:, vcat(5, 7,9, 13,14, 17)])))
    inputTestDB_withDeltaRiTPTN[!, "p(0)"] = pTP_test[:, 1]
    inputTestDB_withDeltaRiTPTN[!, "p(1)"] = pTP_test[:, 2]
    #
## save ##, ouputing trainSet df __ x 19+2 df
savePath = "F:\\UvA\\app\\PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTNandpTP_KNN.csv"
CSV.write(savePath, inputTestDB_withDeltaRiTPTN)


# ==================================================================================================
## count individual sample performance ##
inputTestDB_withDeltaRiTPTN = CSV.read("F:\\UvA\\F\\UvA\\app\\PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTNandpTP_KNN.csv", DataFrame)
sort!(inputTestDB_withDeltaRiTPTN, [:INCHIKEYreal, :INCHIKEY])
    ## define functions for computing weighted P(TP) ##
    function checkID(inID, refID)
        if (inID == refID)
            return true
        else
            return false
        end
    end
    #
    function countID(count, inID, refID)
        acc = count
        if (checkID(inID, refID) == true)
            acc += 1
            return acc
        else
            return 1
        end
    end
    #
    function countP(pro, accPro, inID, refID)
        acc = accPro
        if (checkID(inID, refID) == true)
            acc += pro
            return acc
        else
            return pro
        end
    end
    #
    ## calculate weighted P(TP) for individual sample MS/MS spectrum ##
    count = 1
    colCount = [1]
    accP0 = inputTestDB_withDeltaRiTPTN[1, "p(0)"]
    colP0 = [accP0]
    accP1 = inputTestDB_withDeltaRiTPTN[1, "p(1)"]
    colP1 = [accP1]
    colIDSummary = []
    for id in 1:size(inputTestDB_withDeltaRiTPTN[:, "INCHIKEY"], 1) - 1
        count = countID(count, inputTestDB_withDeltaRiTPTN[id+1, "INCHIKEY"], inputTestDB_withDeltaRiTPTN[id, "INCHIKEY"])
        accP0 = countP(inputTestDB_withDeltaRiTPTN[id+1, "p(0)"], accP0, inputTestDB_withDeltaRiTPTN[id+1, "INCHIKEY"], inputTestDB_withDeltaRiTPTN[id, "INCHIKEY"])
        accP1 = countP(inputTestDB_withDeltaRiTPTN[id+1, "p(1)"], accP1, inputTestDB_withDeltaRiTPTN[id+1, "INCHIKEY"], inputTestDB_withDeltaRiTPTN[id, "INCHIKEY"])
        push!(colCount, count)
        push!(colP0, accP0)
        push!(colP1, accP1)
        if (colCount[id] >= colCount[id+1])
            push!(colIDSummary, 1)
        else
            push!(colIDSummary, 0)
        end
    end
    push!(colIDSummary, 1)
    inputTestDB_withDeltaRiTPTN[!, "countID"] = colCount
    inputTestDB_withDeltaRiTPTN[!, "countP(0)"] = colP0 ./ colCount
    inputTestDB_withDeltaRiTPTN[!, "countP(1)"] = colP1 ./ colCount
    inputTestDB_withDeltaRiTPTN[!, "colIDSummary"] = colIDSummary
    inputTestDB_withDeltaRiTPTN = inputTestDB_withDeltaRiTPTN[inputTestDB_withDeltaRiTPTN."colIDSummary" .== 1, :]

## save ##
savePath = "F:\\UvA\\F\\UvA\\app\\PestMix1-8_test_report_comp_IDs_withDeltaRIandPredictedTPTNandpTP_KNN_ind.csv"
CSV.write(savePath, inputTestDB_withDeltaRiTPTN)
