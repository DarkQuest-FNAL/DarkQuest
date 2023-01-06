#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <TMath.h>

#include "rapidcsv.h"

namespace YieldHelper
{
    const std::string baseDirCSV = "/seaquest/users/yfeng/DarkQuest/DarkQuest/lhe/";

    vector<vector<double>> getBRs(std::string lepName = "Muons")
    {
        // get BR from csv file
        std::string fileName = baseDirCSV + "data/BFto" + lepName + ".txt";
        std::cout << "Reading Branching ratio from " << fileName << std::endl;

        rapidcsv::Document doc(fileName, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams('\t'),
                               rapidcsv::ConverterParams(), rapidcsv::LineReaderParams(true, '#', true));
        std::vector<double> logmasses = doc.GetColumn<double>(0);
        std::vector<double> logbrs = doc.GetColumn<double>(1);

        std::vector<std::vector<double>> brs(logbrs.size(), std::vector<double>(2, 0.));
        for (int i = 0; i < logmasses.size(); i++)
        {
            brs[i][0] = TMath::Power(10, logmasses[i]);
            brs[i][1] = TMath::Power(10, logbrs[i]);
        }

        return brs;
    }

    vector<vector<double>> getDefaultYields(std::string mechName = "Brem")
    {
        // get default yield (i.e., NAp) from csv file
        // mech can be Brem, Eta, or Pion, etc
        // does not depend on the lepton type
        std::string fileName = baseDirCSV + "data/" + mechName + "Yield.txt";
        std::cout << "Reading Default Yields from " << fileName << std::endl;

        rapidcsv::Document doc(fileName, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams('\t'),
                               rapidcsv::ConverterParams(), rapidcsv::LineReaderParams(true, '#', true));
        std::vector<double> logmasses = doc.GetColumn<double>(0);
        std::vector<double> logyields = doc.GetColumn<double>(1);

        std::vector<std::vector<double>> yields(logyields.size(), std::vector<double>(2, 0.));
        for (int i = 0; i < logmasses.size(); i++)
        {
            yields[i][0] = TMath::Power(10, logmasses[i]);
            yields[i][1] = TMath::Power(10, logyields[i]);
        }

        return yields;
    }

    double findValue(std::vector<std::vector<double>> &vec, double mass)
    {
        // find the value in the vector based on the mass
        // linear interpolation if the mass is not in the vector
        if (mass < vec[0][0])
        {
            // std::cout << "mass is smaller than the smallest mass in the vector" << std::endl;
            return vec[0][1];
        }
        else if (mass > vec[vec.size() - 1][0])
        {
            // std::cout << "mass is larger than the largest mass in the vector" << std::endl;
            return vec[vec.size() - 1][1];
        }
        else
        {
            for (int i = 0; i < vec.size() - 1; i++)
            {
                if (mass >= vec[i][0] && mass < vec[i + 1][0])
                {
                    return vec[i][1] + (vec[i + 1][1] - vec[i][1]) / (vec[i + 1][0] - vec[i][0]) * (mass - vec[i][0]);
                }
            }
        }
        return 0.;
    }

    double getYield(std::string mechName, std::string lepName, double mass, double epsilon, bool doPrint = false)
    {
        // get the yield for a given mass
        // mech can be Brem, Eta, or Pion, etc
        // lep can be Muons or Electrons
        std::vector<std::vector<double>> yield0s = getDefaultYields(mechName);
        std::vector<std::vector<double>> brs = getBRs(lepName);
        double yield0 = findValue(yield0s, mass);
        double br = findValue(brs, mass);

        // the yield assumes epsilon = 10^-6, POT = 1.44*10^18,
        // rescale yield as NAp = column2*(epsilon/10^-6)^2 (POT/1.44*10^18)
        const double POT = 1.44e18;
        const double POT0 = 1.44e18;
        const double epsilon0 = 1e-6;

        double yield = yield0 * TMath::Power(epsilon / epsilon0, 2.0) * (POT / POT0);

        if (doPrint)
        {
            std::cout << "for " << mechName << " " << lepName << " mass = " << mass << " epsilon = " << epsilon << std::endl;
            std::cout << " yield x br = " << yield << " x " << br << " = " << yield * br << std::endl;
        }
        return yield * br;
    }
}