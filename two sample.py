
from math import log, sqrt, ceil, pi
import scipy.stats

alpha = float(input("Enter the significance level: "))
assurance = float(input("Enter the assurance probability: "))
theta1 = float(input("Enter the first AUC: "))
theta2 = float(input("Enter the second AUC which is larger than the first AUC: "))
Delta0 = float(input("Enter the lower bound: "))
rho = float(input("Enter the between-AUC correlation: "))
sizeRatio = float(input("Enter the sample size ratio of nondiseased participants to diseased: "))
sdRatio1 = float(input("Enter Test 1's standard deviation ratio of nondiseased participants to diseased: "))
sdRatio2 = float(input("Enter Test 2's standard deviation ratio of nondiseased participants to diseased: "))


beta = 1 - assurance
zalpha = scipy.stats.norm.ppf(1 - alpha / 2)
zbeta = scipy.stats.norm.ppf(1 - beta)

thetaNew = (theta2 - theta1 + 1) / 2

lgtThetaNew = log(thetaNew / (1 - thetaNew))

theta0new = (Delta0 + 1) / 2

lgtTheta0new = log(theta0new / (1 - theta0new))

eta1 = scipy.stats.norm.ppf(theta1) * sqrt(2)
eta2 = scipy.stats.norm.ppf(theta2) * sqrt(2)

fV1 = 0.5 * (scipy.stats.norm.pdf(scipy.stats.norm.ppf(theta1))) ** 2 * (eta1 ** 2 / (2 * (1 + sdRatio1 ** 2) ** 2) * (sizeRatio +
      1 + sdRatio1 ** 4 * (sizeRatio + 1) / sizeRatio) + 2 * (sizeRatio + 1) / (1 + sdRatio1 ** 2) +
      sdRatio1 ** 2 * 2 * (sizeRatio + 1) / (sizeRatio * (1 + sdRatio1 ** 2)))

fV2 = 0.5 * (scipy.stats.norm.pdf(scipy.stats.norm.ppf(theta2))) ** 2 * (eta2 ** 2 / (2 * (1 + sdRatio2 ** 2) ** 2) * (sizeRatio +
      1 + sdRatio2 ** 4 * (sizeRatio + 1) / sizeRatio) + 2 * (sizeRatio + 1) / (1 + sdRatio2 ** 2) +
      sdRatio2 ** 2 * 2 * (sizeRatio + 1) / (sizeRatio * (1 + sdRatio2 ** 2)))

ntotal = ((zbeta + zalpha) / (lgtThetaNew - lgtTheta0new)) ** 2 * 0.25 * (fV1 + fV2 -
       2 * rho * sqrt(fV1) * sqrt(fV2)) / (thetaNew ** 2 * (1 - thetaNew) ** 2)


# make sample sizes integers
nDisease = ceil(pi / 3 * ntotal / (sizeRatio + 1))
nControl = ceil(pi / 3 * sizeRatio * ntotal / (sizeRatio + 1))

ntotal = nDisease + nControl

print("require " + str(ntotal) + " participants with " + str(nDisease) + " diseased + " + str(nControl) + " nondiseased")