//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "MODESWorld.h"

// this is how you setup a parameter in MABE, the function Parameters::register_parameter()takes the
// name of the parameter (catagory-name), default value (which must conform with the type), a the useage message

std::shared_ptr<ParameterLink<int>> MODESWorld::REPORTSPL =
Parameters::register_parameter("WORLD_MODES-REPORTS_terminalReports", 0,
    "should made print reports to the screen? (for debugging) 0: no reports, 1: report only at end of run, 2: report a lot of data every generation");

std::shared_ptr<ParameterLink<int>> MODESWorld::SEXPL =
Parameters::register_parameter("WORLD_MODES-CONTROL_SEX", 0,
    "run with sex? If 0, no sex. If > 0, this value is used for lek size.");

std::shared_ptr<ParameterLink<double>> MODESWorld::TARGETPL =
Parameters::register_parameter("WORLD_MODES-CONTROL_TARGET", 1.0,
    "how much do targets effect fitness?");

std::shared_ptr<ParameterLink<std::string>> MODESWorld::TARGET_LISTSPL =
Parameters::register_parameter("WORLD_MODES-TARGETS_targetLists", std::string("111,101010"),
    "what are the targets of fitness?");

std::shared_ptr<ParameterLink<double>> MODESWorld::FWRPL =
Parameters::register_parameter("WORLD_MODES-CONTROL_FWR", 0.0,
    "how much does being a rare bitsting effect fitness?");

std::shared_ptr<ParameterLink<std::string>> MODESWorld::FWRsizesPL =
Parameters::register_parameter("WORLD_MODES-CONTROL_FWRsizes", std::string("3,6"),
    "how many bit differences are considered when calcluationg fit when rare?\n"
    "if the difference <= first number then the genomes are 'the same'\n"
    "if the difference is >= the second number the genomes are 'different'\n"
    "each bit of observed difference improves score by 1/(max-min)");

std::shared_ptr<ParameterLink<double>> MODESWorld::SHAREPL =
Parameters::register_parameter("WORLD_MODES-CONTROL_SHARE", 0.0,
    "how much does being a rare bitsting effect fitness (using fitness sharing)?");

std::shared_ptr<ParameterLink<double>> MODESWorld::PARASITESPL =
Parameters::register_parameter("WORLD_MODES-CONTROL_PARASITES", 0.0,
    "how much does each parasites effect fitness? if 1 then each parasite drops fitness by selectionBenefit. If 0, no parasites.");

std::shared_ptr<ParameterLink<int>> MODESWorld::popSizePL =
Parameters::register_parameter("WORLD_MODES-BASIC_popSize", 100,
    "population size for main population. If sex, for females. Parasites are seeded one per main pop");


std::shared_ptr<ParameterLink<int>> MODESWorld::generationsPL =
Parameters::register_parameter("WORLD_MODES-BASIC_generations", 20000,
    "number of generation to run model");

std::shared_ptr<ParameterLink<int>> MODESWorld::genomeLengthPL =
Parameters::register_parameter("WORLD_MODES-BASIC_genomeLength", 25,
    "length of main population genomes");

std::shared_ptr<ParameterLink<int>> MODESWorld::genomeLengthSEXPL =
Parameters::register_parameter("WORLD_MODES-SEX_genomeLength", 15,
    "length of sex genomes must be <= main population genomes");

std::shared_ptr<ParameterLink<int>> MODESWorld::genomeLengthPARAPL =
Parameters::register_parameter("WORLD_MODES-PARASITE_genomeLength", 15,
    "length of parasite genomes must be <= main population genomes");

std::shared_ptr<ParameterLink<double>> MODESWorld::mutationRatePL =
Parameters::register_parameter("WORLD_MODES-BASIC_mutationRatePoint", 0.0025,
    "per site chance to flip one bit, only single point mutations are allowed");

std::shared_ptr<ParameterLink<double>> MODESWorld::indelRatePL =
Parameters::register_parameter("WORLD_MODES-BASIC_mutationRateCopyOver", 0.0025,
    "per genome chance to have an copy over mutation (copy than paste that maintains size) of size genome size*[.25,.5]");


std::shared_ptr<ParameterLink<std::string>> MODESWorld::parasiteAltFitnessPL =
Parameters::register_parameter("WORLD_MODES-PARASITE_altFitness", std::string("0"),
    "if not 0 (off), two values (int and double) must be provided.\n"
    "the first value will indicate a number of bits that the parasite must match with the host to allow infection,\n"
    "the second number (double) is the cost for each bit (beyond the first value) that is set to one.");

std::shared_ptr<ParameterLink<double>> MODESWorld::parasiteMigrationRatePL =
Parameters::register_parameter("WORLD_MODES-PARASITE_migrationRate", 0.05,
    "once per host generation each parasite migrates to a new random host with parasiteMigrationRate probablity");

std::shared_ptr<ParameterLink<int>> MODESWorld::parasiteMaxAgePL =
Parameters::register_parameter("WORLD_MODES-PARASITE_maxAge", 20,
    "parasites with age (in parasite generations) > maxAge die... unless they are alone on their host");

std::shared_ptr<ParameterLink<int>> MODESWorld::parasiteGPGPL =
Parameters::register_parameter("WORLD_MODES-PARASITE_gensPerHost", 10,
    "parsites experiance gensPerHost generations (death and birth) for every host generation");

std::shared_ptr<ParameterLink<double>> MODESWorld::parasiteReproCostPL =
Parameters::register_parameter("WORLD_MODES-PARASITE_reproCost", 15.0,
    "cost to produce an offspring");

std::shared_ptr<ParameterLink<int>> MODESWorld::maxParasitesPerHostPL =
Parameters::register_parameter("WORLD_MODES-maxParasitesPerHost", 20,
    "max number of parasites per host. If host has more then extras will be killed of after repro");

std::shared_ptr<ParameterLink<std::string>> MODESWorld::filterDepthListPL =
Parameters::register_parameter("WORLD_MODES-MODES_filterDepthList", std::string("100-,1000/0+,1000/500-"),
    "comma seperated list with clade depth for persistance filter. Each entry will create a seperate MODES tracker.\n"
    "if an entry includes /X then the tracker will start at generation X, default 0\n"
    "+ and - following each value indicate if genomes should be saved");

std::shared_ptr<ParameterLink<double>> MODESWorld::selectionBenefitPL =
Parameters::register_parameter("WORLD_MODES-CONTROL_selectionBenefit", 1.5,
    "selective advantage of each unit increase in score");

std::vector<std::vector<int>> MODESWorld::targets;
double MODESWorld::TARGET;
double MODESWorld::FWR;
double MODESWorld::SHARE;
double MODESWorld::PARASITES;
int MODESWorld::SEX;

MODESWorld::MODESWorld(std::shared_ptr<ParametersTable> PT) : AbstractWorld(PT) {

    // set up targets here so we have access to them for MODES complexity
    std::vector<std::string> tempTargets;
    convertCSVListToVector(TARGET_LISTSPL->get(PT), tempTargets);
    targets.resize(tempTargets.size());
    for (int i = 0; i < tempTargets.size(); i++) {
        for (char c : tempTargets[i]) {
            if (std::isdigit(c)) { // Check if the character is a digit
                // Subtract '0' to convert the character to its integer value
                targets[i].push_back(c - '0');
            }
        }
    }

    TARGET = TARGETPL->get(PT);
    FWR = FWRPL->get(PT);
    SHARE = SHAREPL->get(PT);
    PARASITES = PARASITESPL->get(PT);
    SEX = SEXPL->get(PT);

	popFileColumns.clear();
}

// the evaluate function gets called every generation. evaluate should set values on organisms datamaps
// that will be used by other parts of MABE for things like reproduction and archiving
auto MODESWorld::evaluate(std::map<std::string, std::shared_ptr<Group>>& groups, int analyze, int visualize, int debug) -> void {
    int reports = REPORTSPL->get(PT);

    int lekSize = SEX;

    std::vector<int> FWRsizes;
    convertCSVListToVector(FWRsizesPL->get(PT), FWRsizes);
    if (FWRsizes.size() != 2) {
        std::cout << "  FWRsizes must contain 2 int values. Exitting..." << std::endl;
        exit(1);
    }
    
    double PARASITES = PARASITESPL->get(PT);

    int popSize = popSizePL->get(PT);         // Number of males and females
    int generations = generationsPL->get(PT);    // Run time

    int genomeLength = genomeLengthPL->get(PT);
    int genomeLengthPARA = genomeLengthPARAPL->get(PT);
    int genomeLengthSEX = genomeLengthSEXPL->get(PT);

    double mutationRate = mutationRatePL->get(PT);  // Per org chance to flip 1 bit in either lock or key
    double indelRate = indelRatePL->get(PT);     // Per org chance to add or subtract 1 from the end of lock or key

    std::vector<double> parasiteAltMethod;
    convertCSVListToVector(parasiteAltFitnessPL->get(PT), parasiteAltMethod);
    if (parasiteAltMethod[0] > 0.0) {
        genomeLengthPARA = genomeLength;
    }

    double parasiteMigrationRate = parasiteMigrationRatePL->get(PT);  // Per parasite chance to move to a new host
    int parasiteMaxAge = parasiteMaxAgePL->get(PT);              // Rate per parasite generation at which parasites die
    int parasiteGPG = parasiteGPGPL->get(PT);                 // Number of parasite generations per host generation
    double parasiteReproCost = parasiteReproCostPL->get(PT);           // Cost for a parasite to make a baby
    int maxParasitesPerHost = maxParasitesPerHostPL->get(PT);

    double selectionBenefit = selectionBenefitPL->get(PT);

    // add some MODES trackers!
    std::vector<std::string> tempTrackerData;
    convertCSVListToVector(filterDepthListPL->get(PT), tempTrackerData);
    std::vector<MODES_TRACKER> trackers;
    int filterMax = 0;
    for (auto& elem : tempTrackerData) {
        std::vector<std::string> trakerParts;
        bool saveGenomes = (elem.back() == '+') ? true : false;
        elem.pop_back(); // remove char that indicates saving genome
        convertCSVListToVector(elem, trakerParts,'/');
        int num = std::stoi(trakerParts[0]);
        int trackerStart;
        trackerStart = (trakerParts.size() == 1) ? 0 : stoi(trakerParts[1]);
        trackers.emplace_back(std::stoi(elem), saveGenomes, popSize, trackerStart);
        filterMax = std::max(filterMax, std::stoi(elem));
    }

    //int changeDelay = changeDelayPL->get(PT);

    // Constructing the EXPERIMENT and DATANAME strings
    std::ostringstream experiment_stream;
    experiment_stream << "SEX_" << SEX << "_FWR_" << FWR << "_PAR_" << PARASITES
        << "_PS_" << popSize << "_mr_" << mutationRate << "_ir_" << indelRate;
    std::string EXPERIMENT = experiment_stream.str();

    std::string DATANAME = EXPERIMENT + "_data/";


    // Initialize populations
    std::cout << "Initializing populations of size " << popSize << "\n";
    std::vector<std::vector<int>> mPop;
    mPop.reserve(popSize); // Reserve space for popSize elements for efficiency

    // Vectors for phylogeny data
    std::vector<int> phylo_m_current(popSize); // used to track phylogony
    //std::iota(phylo_m_current.begin(), phylo_m_current.end(), 0); // Fill with 0, 1, 2, ...
    std::vector<int> phylo_m_next(popSize);    // used to track phylogony

    for (int i = 0; i < popSize; ++i) {
        // Create a bitstring (vector<int>) of this length, initializing it with random values
        if (i < 10) {
            std::vector<int> bitString(genomeLength);
            std::generate(bitString.begin(), bitString.end(), []() { return Random::getInt(0, 1); });

            // Add the bitstring to mPop
            mPop.push_back(bitString);
            phylo_m_current[i] = i;
        }
        else { // this will start us out with 10 random genomes. This avoids too many types for FWR
            int pick = Random::getIndex(10);
            mPop.push_back(mPop[pick]);
            phylo_m_current[i] = pick;
        }
    }
    std::vector<std::vector<int>> next_mPop(popSize);


    std::vector<std::vector<int>> mPop_locks;
    mPop_locks.reserve(popSize);
    for (int i = 0; i < popSize; ++i) {
        if (i < 10) {
            std::vector<int> bitString(genomeLengthSEX);
            std::generate(bitString.begin(), bitString.end(), []() { return Random::getInt(0, 1); });
            mPop_locks.push_back(bitString);
        }
        else { // this will start us out with 10 random genomes. This avoids too many types for FWR
            int pick = Random::getIndex(10);
            mPop_locks.push_back(mPop_locks[pick]);
        }
    }
    std::vector<std::vector<int>> next_mPop_locks(popSize);


    //std::vector<std::vector<int>> fPop;
    //fPop.reserve(popSize); // Reserve space for popSize elements for efficiency
    //
    //for (int i = 0; i < popSize; ++i) {
    //    std::vector<int> bitString(Random::getInt(initGenomeLengthSEX[0], initGenomeLengthSEX[1]));
    //    std::generate(bitString.begin(), bitString.end(), []() { return Random::getInt(0, 1); });
    //    fPop.push_back(bitString);
    //}
    //std::vector<std::vector<int>> next_fPop(popSize);
    
    std::vector<std::vector<Parasite>> pPop(popSize);
    int i = 0;
    for (auto& vec : pPop) {
        vec.reserve(1); // Reserve space for 1 Parasite
        std::vector<int> genome(genomeLengthPARA);
        if (i < 10) {
            if (parasiteAltMethod[0] == 0) {
                std::generate(genome.begin(), genome.end(), []() { return Random::getInt(0, 1); });
            }
            else { // is using alt method, start all parasites as all 0
                std::generate(genome.begin(), genome.end(), []() { return 1; });
            }
            i++;
        }
        else {
            int pick = Random::getIndex(10);
            genome = pPop[pick][0].genome;
        }
        vec.emplace_back(std::move(genome), Random::getDouble(0, parasiteReproCost), Random::getInt(0, parasiteMaxAge));
    }

    std::vector<std::vector<Parasite>> next_pPop(popSize);

    // Vectors to hold counts of unique types for each generation
    std::vector<int> m_types_counts;
    std::vector<int> f_types_counts;
    std::vector<int> p_types_counts;

    // counts of total number of parasites per generation
    std::vector<int> p_counts;

    //std::vector<std::vector<int>> m_lengths_lists;  // list of list of mPop genome lengths per generation
    //std::vector<std::vector<int>> f_lengths_lists;  // list of list of fPop genome lengths + mPop_locks lengths per generation
    //std::vector<std::vector<int>> p_lengths_lists;  // list of list of parasites genome lengths per generation
    //std::vector<int> max_m_lengths; // list of longest genome in mPop every generation
    //std::vector<int> max_f_lengths; // list of longest genome in fPop and mPop_locks every generation
    //std::vector<int> max_p_lengths; // list of longest genome in pPop every generation

    std::vector<double> scores(popSize);

    // Parameters for file management
    int generationsPerFile = 10000; // Number of generations per file
    if (generationsPerFile > generations) {
        generationsPerFile = generations;
    }
    bool saveRunTimeStats = true;

    std::string filename_genomes = "BLANK";
    std::string filename_MODES = "MODES_data.csv";
    std::string filename_runTimeData = "runTimeStats.csv";
    std::string runTimeData_headerString = "BLANK";

    int currentFileStartGen = 0; // Starting generation of the current file

    //std::vector<std::vector<int>> mPopHistory(popSize);// mPop);
    //std::vector<int> phylo_m_history(phylo_m_current);
    //std::vector<double> score_m_history(popSize);

    std::cout << "starting up MODES World: " << std::endl;
    std::cout << "   FWR................... " << FWR << std::endl;
    std::cout << "   SHARE................. " << SHARE << std::endl;
    std::cout << "     FWR/SHARE sizes..... " << FWRsizes[0] << ", " << FWRsizes[1] << std::endl;
    std::cout << "   SEX................... " << SEX << std::endl;
    std::cout << "   PARASITES............. " << PARASITES << std::endl;
    std::cout << "   TARGET................ " << TARGET << std::endl;
    std::cout << "\n   global states:" << std::endl;
    std::cout << "      popSize............ " << popSize << std::endl;
    std::cout << "      generations........ " << generations << std::endl;
    std::cout << "      selectionBenefit... " << selectionBenefit << std::endl;
    std::cout << "\n   genomes:" << std::endl;
    std::cout << "      genomeLength....... " << genomeLength << std::endl;
    std::cout << "      point rate......... " << mutationRate << std::endl;
    std::cout << "      indel rate......... " << indelRate << std::endl;
    std::cout << "\n   parasites:" << std::endl;
    std::cout << "      parasiteAltMethod.. " << ((parasiteAltMethod[0] == 0) ? std::string("OFF") : std::to_string(parasiteAltMethod[0]) + "/" + std::to_string(parasiteAltMethod[1])) << std::endl;
    std::cout << "      genomeLength....... " << genomeLengthPARA << ((parasiteAltMethod[0] == 0) ? std::string("") : std::string("  override by alt method")) << std::endl;
    std::cout << "      MigrationRate...... " << parasiteMigrationRate << std::endl;
    std::cout << "      MaxAge............. " << parasiteMaxAge << std::endl;
    std::cout << "      PG/HG.............. " << parasiteGPG << std::endl;
    std::cout << "      reproCost.......... " << parasiteReproCost << std::endl;
    std::cout << "      maxParasitesPerHost " << maxParasitesPerHost << std::endl;
    std::cout << "\n   sexual selection:" << std::endl;
    std::cout << "      genomeLength....... " << genomeLengthSEX << std::endl;
    std::cout << "      lekSize............ " << lekSize << std::endl;
    std::cout << "\n   data:" <<  std::endl;
    std::cout << "      report level....... " << reports << std::endl;
    std::cout << std::endl;
    std::cout << "   found targets : " << std::endl;
    for (auto t : targets) {
        std::cout << "      ";
        for (auto v : t) {
            std::cout << v;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "   MODES:" << std::endl;
    int tr_counter = 0;
    for (auto& tr : trackers) {
        std::cout << "      " << tr.name << ":" << std::endl;
        std::cout << "        filterDepth........ " << tr.filterDepth << std::endl;
        std::cout << "        change delay....... " << tr.changeDelay << std::endl;
        std::cout << "        save genomes....... " << tr.saveGenomes << std::endl;
    }
    std::cout << std::endl;

    for (int gen = 0; gen < generations + 1 + filterMax; ++gen) {

        // report run status to terminal
        if (gen % 100 == 0) {
            if (reports > 1) {
                std::cout << "population:" << std::endl;
                for (auto m : mPop) {
                    for (auto b : m)  std::cout << b;
                    std::cout << std::endl;
                }
                if (SEX) {
                    std::cout << "females:" << std::endl;
                    for (auto f : mPop_locks) {
                        for (auto b : f) std::cout << b;
                        std::cout << std::endl;
                    }
                }
            }
            std::string pc = (p_counts.size()>0)?std::to_string(p_counts.back()):std::string("-");
            std::cout << "generation: " << gen << "  parasitesCount: " << pc << std::endl;
        }


        //////////////////////////////////////
        // start data tracking
        //////////////////////////////////////

        if (gen % generationsPerFile == 0) {

            // save run time stats
            if (saveRunTimeStats) {
                std::string RUNTIMEDATA = "";
                if (gen == 0) { // if it's gen 0, setup header string
                    if (SEX && PARASITES) {
                        //runTimeData_headerString = "m_types,max_m_length,m_lengths_list,f_types,max_f_length,f_lengths_list,p_types,max_p_length,p_lengths_list,p_count";
                        runTimeData_headerString = "m_types,f_types,p_types,p_count";
                    }
                    if (SEX && !PARASITES) {
                        //runTimeData_headerString = "m_types,max_m_length,m_lengths_list,f_types,max_f_length,f_lengths_list";
                        runTimeData_headerString = "m_types,f_types";
                    }
                    if (!SEX && PARASITES) {
                        //runTimeData_headerString = "m_types,max_m_length,m_lengths_list,p_types,max_p_length,p_lengths_list,p_count";
                        runTimeData_headerString = "m_types,p_types,p_count";
                    }
                    if (!SEX && !PARASITES) {
                        //runTimeData_headerString = "m_types,max_m_length,m_lengths_list";
                        runTimeData_headerString = "m_types";
                    }
                }
                else { // not gen 0, open file in "append" mode and add some data
                    for (int i = 0; i < m_types_counts.size(); i++) {
                        if (SEX && PARASITES) {
                            RUNTIMEDATA += std::to_string(m_types_counts[i]) + ",";
                            RUNTIMEDATA += std::to_string(f_types_counts[i]) + ",";
                            RUNTIMEDATA += std::to_string(p_types_counts[i]) + ",";
                            RUNTIMEDATA += std::to_string(p_counts[i]) + "\n";
                        }

                        if (SEX && !PARASITES) {
                            RUNTIMEDATA += std::to_string(m_types_counts[i]) + ",";
                            RUNTIMEDATA += std::to_string(f_types_counts[i]) + "\n";
                        }
                        if (!SEX && PARASITES) {
                            RUNTIMEDATA += std::to_string(m_types_counts[i]) + ",";
                            RUNTIMEDATA += std::to_string(p_types_counts[i]) + ",";
                            RUNTIMEDATA += std::to_string(p_counts[i]) + "\n";
                        }
                        if (!SEX && !PARASITES) {
                            RUNTIMEDATA += std::to_string(m_types_counts[i]) + "\n";
                        }

                        // old version with sizes
                        /*
                        if (0) {
                            if (SEX && PARASITES) {
                                RUNTIMEDATA += std::to_string(m_types_counts[i]) + "," + std::to_string(max_m_lengths[i]) + ",\"";
                                for (auto v : m_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\",";

                                RUNTIMEDATA += std::to_string(f_types_counts[i]) + "," + std::to_string(max_f_lengths[i]) + ",\"";
                                for (auto v : f_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\",";

                                RUNTIMEDATA += std::to_string(p_types_counts[i]) + "," + std::to_string(max_p_lengths[i]) + ",\"";
                                for (auto v : p_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\",";
                                RUNTIMEDATA += std::to_string(p_counts[i]) + "\n";
                            }

                            if (SEX && !PARASITES) {
                                RUNTIMEDATA += std::to_string(m_types_counts[i]) + "," + std::to_string(max_m_lengths[i]) + ",\"";
                                for (auto v : m_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\",";

                                RUNTIMEDATA += std::to_string(f_types_counts[i]) + "," + std::to_string(max_f_lengths[i]) + ",\"";
                                for (auto v : f_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\"\n";
                            }
                            if (!SEX && PARASITES) {
                                RUNTIMEDATA += std::to_string(m_types_counts[i]) + "," + std::to_string(max_m_lengths[i]) + ",\"";
                                for (auto v : m_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\",";

                                RUNTIMEDATA += std::to_string(p_types_counts[i]) + "," + std::to_string(max_p_lengths[i]) + ",\"";
                                for (auto v : p_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\",";
                                RUNTIMEDATA += std::to_string(p_counts[i]) + "\n";
                            }
                            if (!SEX && !PARASITES) {
                                RUNTIMEDATA += std::to_string(m_types_counts[i]) + "," + std::to_string(max_m_lengths[i]) + ",\"";
                                for (auto v : m_lengths_lists[i]) RUNTIMEDATA += std::to_string(v) + ",";
                                RUNTIMEDATA += "\"\n";
                            }
                        } // old vserion with sizes
                        */


                    }
                }

                FileManager::openAndWriteToFile(filename_runTimeData, RUNTIMEDATA, runTimeData_headerString);
                FileManager::closeFile(filename_runTimeData);

                // clear out the stats
                m_types_counts.clear();
                f_types_counts.clear();
                p_types_counts.clear();

                p_counts.clear();

                //m_lengths_lists.clear();
                //f_lengths_lists.clear();
                //p_lengths_lists.clear();
                //max_m_lengths.clear();
                //max_f_lengths.clear();
                //max_p_lengths.clear();

            } // done save runtime stats

        }

        //////////////////////////////////////
        // done data tracking
        //////////////////////////////////////
        

        // get unique type counts and counts by length for main population
        //std::vector<int> m_lengths(2 * maxGenomeLength + 1, 0);
        std::map<std::vector<int>, int> m_count_list;
        for (const auto& m : mPop) {
            //m_lengths[std::min(2 * maxGenomeLength, (int)m.size())]++;
            m_count_list[m]++;
        }
        //m_lengths_lists.push_back(m_lengths);       // Push counts by length for this generation
        m_types_counts.push_back(m_count_list.size());   // Push back the number of unique genomes

        // Find the maximum length for this generation
        //int max_length = 0;
        //for (int i = m_lengths.size() - 1; i >= 0; --i) {
        //    if (m_lengths[i] > 0) {
        //        max_length = i;
        //        break;
        //    }
        //}
        //max_m_lengths.push_back(max_length);

        if (SEX) { // sex setup and data collection
            //auto combinedLocks = fPop;
            //combinedLocks.insert(combinedLocks.end(), mPop_locks.begin(), mPop_locks.end());

            // get unique type counts and counts by length
            //std::vector<int> f_lengths(2 * maxGenomeLength + 1, 0);
            std::map<std::vector<int>, int> f_count_list;
            for (const auto& f : mPop_locks) {
            //    f_lengths[std::min(2 * maxGenomeLength,(int)f.size())]++;
                f_count_list[f]++;
            }
            //f_lengths_lists.push_back(f_lengths);       // Push counts by length for this generation
            f_types_counts.push_back(f_count_list.size());   // Push back the number of unique genomes

            // Find the maximum length for this generation
            //int max_length = 0;
            //for (int i = f_lengths.size() - 1; i >= 0; --i) {
            //    if (f_lengths[i] > 0) {
            //        max_length = i;
            //        break;
            //    }
            //}
            //max_f_lengths.push_back(max_length);

        } // end SEX setup

        scores = getFit(mPop, SHARE, FWR, FWRsizes, TARGET, targets);

        if (PARASITES) {
            //// get unique type counts and counts by length for parasites
            //std::vector<int> p_lengths(2 * maxGenomeLength + 1, 0);
            std::map<std::vector<int>, int> p_count_list;
            int p_count = 0;
            for (const auto& host : pPop) {
                for (const auto& parasite : host) {
            //        p_lengths[std::min(2 * maxGenomeLength,(int)parasite.genome.size())]++;
                    p_count_list[parasite.genome]++;
                    p_count++;
                }
            }
            //p_lengths_lists.push_back(p_lengths);       // Push counts by length for this generation
            p_types_counts.push_back(p_count_list.size());   // Push back the number of unique genomes
            
            // Find the maximum length for this generation
            //int max_length = 0;
            //for (int i = p_lengths.size() - 1; i >= 0; --i) {
            //    if (p_lengths[i] > 0) {
            //        max_length = i;
            //        break;
            //    }
            //}
            //max_p_lengths.push_back(max_length);

            p_counts.push_back(p_count);

            // migrate parasites
            std::vector<Parasite>  interHost;     //put migrating parasites here so they don't get moved more than once
            for (auto& host : pPop) { //for each host, migrate some parasites to interHosts
                if ((int)host.size() > 1) {  // a parasite will not migrate if they are alone - this is to avoid extiction, also they are happy when no one is around...
                    for (int r = Random::getBinomial(((int)host.size()) - 1, parasiteMigrationRate); r > 0; r--) { // how many parasites move from this host?
                        int whichParasite = Random::getIndex(host.size());
                        interHost.push_back(host[whichParasite]);
                        host[whichParasite] = host.back();
                        host.pop_back();
                    }
                }
            } // end migrate parasites

            // migrant parasites need new homes
            for (auto parasite : interHost) {
                pPop[Random::getIndex(popSize)].push_back(parasite); // for each parasite that's migrating, put them on a random host
            }

            // evolve parasites
            int host_index = -1;
            for (auto& host : pPop) {
                host_index++;
                
                bool pReport = gen % 100 == 0 && host_index == 0;
                pReport = false; // comment out to show a lot of parasite related runtime data

                for (int parasite_index = 0; parasite_index < host.size(); parasite_index++) { // clear out resource
                    host[parasite_index].resource = 0;
                    host[parasite_index].age = 0;
                }
                double hostLosses = 0.0;
                for (int p_gen = 0; p_gen < parasiteGPG; p_gen++) {
                    // kill parasites from old age
                    for (int parasite_index = host.size() - 1; parasite_index >= 0; parasite_index--) {
                        if (host[parasite_index].age > parasiteMaxAge && host.size() > 1 && Random::P(.5)) {
                            host[parasite_index] = host.back();
                            host.pop_back();
                        }
                        else {
                            host[parasite_index].age++;
                        }
                    }

                    // next get score for each parasite on this host.the better the bitstring match, the higher the score.
                    // parasites total score can not exceed maxTheft, if it does, scores will be reduced
                    std::vector<double> pScores(host.size(), 0.0);
                    for (int parasite_index = 0; parasite_index < host.size(); parasite_index++) {
                        if (parasiteAltMethod[0] == 0) {
                            double hLen = (double)mPop[host_index].size();
                            double pLen = (double)host[parasite_index].genome.size();

                            //////////////////////////////////////////////
                            //////////////////////////////////////////////
                            //////////////////////////////////////////////
                            //////////////////////////////////////////////

                            //pScores[parasite_index] = std::pow((10.0 - static_cast<double>(findMinMismatches(mPop[host_index], host[parasite_index].genome))) / 10.00, 2.0);
                            auto [matchScore, temp] = findBestMatch(encodeGenome(mPop[host_index]), encodeGenome(host[parasite_index].genome));
                            //double halfParasiteGenome = static_cast<double>(genomeLengthPARA) / 1.0;
                            //pScores[parasite_index] = std::pow(std::max(0.0, static_cast<double>(matchScore) - 0) / halfParasiteGenome, 1.0);
                            pScores[parasite_index] = std::pow(static_cast<double>(matchScore) / static_cast<double>(genomeLengthPARA), 2.0);
                            //std::cout << "genome:  " << u642s(encodeGenome(mPop[host_index])) << "\ntarget:  " << u642s(encodeGenome(host[parasite_index].genome)) << "\nmatch:   " << u642s(temp) << "\n         " << std::endl << std::endl;;
                            //std::cout << pScores[parasite_index] << std::endl;
                        }
                        else {
                            auto encoded_host = encodeGenome(mPop[host_index]);
                            if (__builtin_popcountll(encoded_host) < 1 + parasiteAltMethod[0]) { // if host trys to escape "down" parasite wins!
                                pScores[parasite_index] = 1.0;
                            }
                            else {
                                auto encoded_parasite = encodeGenome(host[parasite_index].genome);
                                pScores[parasite_index] = std::max(.01,
                                    // note: -1 and -1.0 are to deal with the leading 1 in the genome encoding
                                    static_cast<double>(__builtin_popcountll(encoded_host & encoded_parasite) - 1 >= parasiteAltMethod[0]) - // do they match enough bits
                                    (std::max(0.0, static_cast<double>(__builtin_popcountll(encoded_parasite)) - (parasiteAltMethod[0] + 1.0)) * parasiteAltMethod[1])
                                );         // lose cost for each extra, beyond the min
                                //std::cout << "genome:  " << u642s(encodeGenome(mPop[host_index])) << "\ntarget:  " << u642s(encodeGenome(host[parasite_index].genome)) << "\n         " << u642s(encodeGenome(mPop[host_index]) & encoded_parasite) << "\n         " << std::endl << std::endl;;
                            }
                        //std::cout << __builtin_popcountll(encodeGenome(mPop[host_index]) & encoded_parasite) << " " << parasiteAltMethod[0] << " " << static_cast<double>(__builtin_popcountll(encodeGenome(mPop[host_index]) & encoded_parasite) >= parasiteAltMethod[0]) << "  score: " << pScores[parasite_index] << std::endl;
                        }

                        //////////////////////////////////////////////
                        //////////////////////////////////////////////
                        //////////////////////////////////////////////
                        //////////////////////////////////////////////

                    }

                    hostLosses += host.size();

                    // distribute resources - parasites resouce increase by their pScore
                    for (int parasite_index = 0; parasite_index < host.size(); parasite_index++) {
                        if (pReport) std::cout << "      score:" << pScores[parasite_index] << " r:" << host[parasite_index].resource;
                        host[parasite_index].resource += Random::getDouble(.9,1.1)*pScores[parasite_index];
                        if (pReport) std::cout << " >>> " << host[parasite_index].resource << std::endl;
                    }
                    
                    // now birth some parasites - if parasite resource[1] is > parasiteReproCost make a baby
                    int temp_num_parasites = host.size();
                    for (int parasite_index = 0; parasite_index < temp_num_parasites; parasite_index++) {
                        if (host[parasite_index].resource > parasiteReproCost) {
                            host.push_back(Parasite( mutate(host[parasite_index].genome, mutationRate, indelRate, 3), 0.0, 0 )); // the wee ones start with no resource, and age 0
                            host[parasite_index].resource -= parasiteReproCost;  // kids cost resource
                        }
                    }

                    // kill parasites from overcrowding
                    while (host.size() > maxParasitesPerHost) {
                        host[Random::getIndex(host.size())] = host.back();
                        host.pop_back();
                    }

                } //  end parasite gens (p_gens)

                // this host is done, now we just need to steal the resources
                if (pReport) std::cout << gen << " para_count:" << host.size() << "    hostLosses:" << hostLosses << " parasiteGPG:" << parasiteGPG << " scores[host_index]:" << scores[host_index];
                scores[host_index] -= (PARASITES * (hostLosses / (double)parasiteGPG));
                if (pReport) std::cout << " >>> " << scores[host_index] << std::endl;

            } // end this host
        } // end PARASITES

        double maxScore = *std::max_element(scores.begin(), scores.end());
        //std::cout << "scores: ";
        for (auto& score : scores) {
            //std::cout << score << " / ";
            score = std::pow(selectionBenefit, score - maxScore);
            //std::cout << score << "  ";
        }
        //std::cout << std::endl;

        for (auto& tr : trackers) {
            tr.updateMODES(mPop, gen, mPop_locks);
            tr.updateScoreHistory(scores, gen);
        }

        int newPopIndex = 0;

        while (newPopIndex < popSize) {

            std::vector<double> cumulative_scores;

            if (SEX) {
                int f_index = Random::getIndex(popSize); // select a female at random
                
                std::vector<int> lek = roulette_wheel_selection(scores, cumulative_scores, lekSize);

                int m_index = -1;
                double bestMatch = 0;

                for (int i = 0; i < lekSize; i++) {
                    auto [matchScore,temp] = findBestMatch(encodeGenome(mPop[lek[i]]), encodeGenome(mPop_locks[f_index]));
                    if (matchScore > bestMatch) {
                        m_index = lek[i];
                        bestMatch = matchScore;
                    }
                }

                if (m_index != -1) { // if a male was selected
                    auto whichLock = Random::getInt(0, 1) ? mPop_locks[m_index] : mPop_locks[f_index]; // does the lock come from mom or dad?
                    //next_fPop[newPopIndex] = mutate(whichLock, mutationRate, indelRate, sizeRate, 3, maxGenomeLengthSEX);
                    next_mPop_locks[newPopIndex] = mutate(whichLock, mutationRate, indelRate, 3);

                    next_mPop[newPopIndex] = mutate(mPop[m_index], mutationRate, indelRate, 3);
                    

                    for (auto& tr : trackers) {
                        tr.updateNextPhylo(m_index, newPopIndex);
                    }

                    if (PARASITES) {
                        next_pPop[newPopIndex] = pPop[m_index];
                    }
                    newPopIndex++;
                }
            }
            else{ // not SEX
                auto picks = roulette_wheel_selection(scores, cumulative_scores);;
                for (auto pick : picks) {
                    next_mPop[newPopIndex] = mutate(mPop[pick], mutationRate, indelRate, 3);
                    

                    for (auto& tr : trackers) {
                        tr.updateNextPhylo(pick, newPopIndex);
                    }

                    if (PARASITES) {
                        next_pPop[newPopIndex] = pPop[pick];
                    }
                    newPopIndex++;
                }
            }
        }


        for (auto& tr : trackers) {
            tr.updatePhylo();
        }

        mPop = next_mPop;

        if (PARASITES) {
            pPop = next_pPop;
        }

        if (SEX) {
            //fPop = next_fPop;
            mPop_locks = next_mPop_locks;
        }

    } // end generations

    std::cout << "\n\n-- done --" << std::endl;

    if (reports == 1) {
        std::cout << "FINAL REPORT" << std::endl;
        std::cout << "population:" << std::endl;
        for (auto m : mPop) {
            for (auto b : m) std::cout << b;
            std::cout << std::endl;
        }
        if (SEX) {
            std::cout << "females:" << std::endl;
            for (auto f : mPop_locks) {
                for (auto b : f) std::cout << b;
                std::cout << std::endl;
            }
        }
        if (PARASITES) {
            std::cout << "Parasite count: " << p_counts.back() << std::endl;
        }
    }
    
    exit(0); // everything ran correctly!
} // end evaluate

// the requiredGroups function lets MABE know how to set up populations of organisms that this world needs
auto MODESWorld::requiredGroups() -> std::unordered_map<std::string, std::unordered_set<std::string>> {
	return { { "root::", { } } };
}
