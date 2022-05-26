
from math import log, sqrt, ceil, pi
import scipy.stats

alpha = float(input("Enter the significance level: "))
assurance = float(input("Enter the assurance probability: "))
theta = float(input("Enter the AUC: "))
theta0 = float(input("Enter the lower bound: "))
sizeRatio = float(input("Enter the sample size ratio of nondiseased participants to diseased: "))
sdRatio = float(input("Enter the standard deviation ratio of nondiseased participants to diseased: "))

beta = 1 - assurance

zalpha = scipy.stats.norm.ppf(1 - alpha / 2)
zbeta = scipy.stats.norm.ppf(1 - beta)

lgtTheta = log(theta / (1 - theta))
lgtTheta0 = log(theta0 / (1 - theta0))

eta = scipy.stats.norm.ppf(theta) * sqrt(2)

fV = 0.5 * (scipy.stats.norm.pdf(scipy.stats.norm.ppf(theta))) ** 2 * (eta ** 2 / (2 * (1 +
        sdRatio ** 2) ** 2) * (sizeRatio + 1 + sdRatio ** 4 * (sizeRatio + 1) / sizeRatio) +
        2 * (sizeRatio + 1) / (1 + sdRatio ** 2) +
        sdRatio ** 2 * 2 * (sizeRatio + 1) / (sizeRatio * (1 + sdRatio ** 2)))

ntotal = ((zbeta + zalpha) / (lgtTheta - lgtTheta0)) ** 2 * fV / (theta ** 2 * (1 - theta) ** 2)

# make sample size integers
nDisease = ceil(pi / 3 * ntotal / (sizeRatio + 1))
nControl = ceil(pi / 3 * sizeRatio * ntotal / (sizeRatio + 1))
ntotal = nDisease + nControl

print("require " + str(ntotal) + " participants with " + str(nDisease) + " diseased + " + str(nControl) + " nondiseased")