//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once    // directive to insure that this .h file is only included one time

#include <World/AbstractWorld.h> // AbstractWorld defines all the basic function templates for worlds
#include <string>
#include <memory> // shared_ptr
#include <map>

#include <limits.h>

#include <vector>
#include <algorithm>
#include <functional>  // for std::hash and std::function
#include <numeric>     // for std::accumulate
#include <cmath>  // for pow, log2 and fabs
#include <bitset>
#include <queue>


//encodeGenome(originalGenome);
//decodeGenome(encoded);
//std::cout << "Encoded (as bits): " << std::bitset<64>(encoded) << std::endl;
//std::unordered_set<uint64_t> fridge; // longterm storage for novelty metric
//std::queue<std::unordered_set<uint64_t>> fridgeQueue; // stores genomes that have been collected but are not yet in the fridge


class MODESWorld : public AbstractWorld {

public:

    static std::shared_ptr < ParameterLink<int>> REPORTSPL;
    //static std::shared_ptr < ParameterLink<bool>> saveGenomesPL;
    static std::shared_ptr < ParameterLink<double>> TARGETPL;
    static std::shared_ptr < ParameterLink<std::string>> TARGET_LISTSPL;
    static std::shared_ptr < ParameterLink<int>> SEXPL;
    static std::shared_ptr < ParameterLink<double>> FWRPL;
    static std::shared_ptr < ParameterLink<double>> SHAREPL;
    static std::shared_ptr < ParameterLink<std::string>> FWRsizesPL;
    static std::shared_ptr < ParameterLink<double>> PARASITESPL;

    static std::shared_ptr < ParameterLink<int>> popSizePL;
    static std::shared_ptr < ParameterLink<int>> generationsPL;

    static std::shared_ptr < ParameterLink<int>> genomeLengthPL;

    static std::shared_ptr < ParameterLink<int>> genomeLengthSEXPL;

    static std::shared_ptr < ParameterLink<int>> genomeLengthPARAPL;

    static std::shared_ptr < ParameterLink<double>> mutationRatePL;
    static std::shared_ptr < ParameterLink<double>> indelRatePL;

    static std::shared_ptr < ParameterLink<std::string>> parasiteAltFitnessPL;
    static std::shared_ptr < ParameterLink<double>> parasiteMigrationRatePL;
    static std::shared_ptr < ParameterLink<int>> parasiteMaxAgePL; 
    static std::shared_ptr < ParameterLink<int>> parasiteGPGPL; 
    static std::shared_ptr < ParameterLink<double>> parasiteReproCostPL;
    static std::shared_ptr < ParameterLink<int>> maxParasitesPerHostPL;

    static std::shared_ptr < ParameterLink<std::string>> filterDepthListPL;
//    static std::shared_ptr < ParameterLink<int>> changeDelayPL;

    static std::shared_ptr < ParameterLink<double>> selectionBenefitPL;


    static std::vector<std::vector<int>> targets;
    static double TARGET;
    static double FWR;
    static double SHARE;
    static double PARASITES;
    static int SEX;

    // IVY: returns the length of a uint64_t genome

    // IVY: Generate a mask with offset and size (ex. 1110 has offset 3 and size 3)
    // perhaps should add some error checking later?
    static uint64_t generateMask(int offset, int size) {
        uint64_t copyMask = (1 << (offset + 1)) - 1;

        // Now need to chop off end of mask?
        if (offset + 1 > size)
        {
            copyMask ^= (1 << (offset - size + 1)) - 1;
        }
        return copyMask;
    }

    static std::string g2s(std::vector<int> genome) {
        std::string s = "";
        for (auto& v : genome) s += std::to_string(v);
        return(s);
    }

    static std::string u642s(uint64_t bits) {
        std::string result;
        for (int i = 63; i >= 0; --i) { // Start from the most significant bit
            result += ((bits >> i) & 1) ? '1' : '0';
        }
        return result;
    }

    // Encodes a genome sequence into an uint64_t, assuming an implicit leading 1
    static uint64_t encodeGenome(const std::vector<int>& genome) {
        uint64_t encoded = 1ULL; // Start with the marker bit set
        for (int bit : genome) {
            encoded = (encoded << 1) | bit;
        }
        return encoded;
    }

    // Decodes an encoded genome back into a vector<int>, skipping the marker bit
    static std::vector<int> decodeGenome(uint64_t encoded) {
        std::vector<int> genome;
        bool start = false;
        for (int i = 63; i >= 0; --i) {
            int bit = (encoded >> i) & 1;
            if (bit == 1 && !start) { // Find the marker bit
                start = true; // Start recording bits after this
                continue;
            }
            if (start) {
                genome.push_back(bit);
            }
        }
        return genome;
    }

    /**
     * matchBits calculates the maximum similarity ratio between a key and a circular lock.
     * This function compares a key vector to a lock vector to find the highest ratio of matching bits to the key's length.
     * The lock is treated as circular, allowing the key to be aligned starting at any position in the lock.
     *
     * @param key A vector of integers representing the key.
     * @param lock A vector of integers representing the circular lock.
     * @return A double representing the highest similarity ratio found, ranging from 0.0 (no match) to 1.0 (perfect match).
     * If the key is longer than the lock or if either the key or lock is empty, the function returns 0.0.
     */

    // Implementation using uint64_t instead of vectors
    static double matchBits(const uint64_t key, uint64_t lock) {
        // Make sure they're both nonempty
        // This is important to using __builtin_clz later (behavior undefined for 0)
        // Even if "empty", they'll still have a 1 bit set to signify start
        if (key == 1 || lock == 1) {
            return 0.0;
        }

        // Check whether key has less leading zeros (aka key is bigger than lock)
        // Not super suer whether the ll is needed, but adding for safety

        int keySize = calculateDataLength(key); // The additional -1 is for the leading 1
        int lockSize = calculateDataLength(lock);

        if (keySize > lockSize) {
            return 0.0;
        }

        int bestMatch = 0;

        // IVY:
        // TODO: There is hacky way to add difference onto front of lock and then iterate through, but
        // I'm not sure if we'll have the hard cap on genome length

        // Rough idea is to take xor then count number of bits with __builtin_popcountll()
    }

    static double matchBits(const std::vector<int>& key, const std::vector<int>& lock) {
        if (key.empty() || lock.empty()) {  // Validate non-empty input vectors
            return 0.0;
        }
        if (key.size() > lock.size()) {  // lock is too small for key
            return 0.0;
        }

        int bestMatch = 0;
        // Loop around the lock to compare the key against all possible starting positions
        for (size_t offset = 0; offset < lock.size(); ++offset) {  // Iterate over the lock to find the best match for the key
            // Count the number of matching bits between the key and the corresponding section of the lock
            int current_match = 0;
            // Count matching bits between the key and the lock at the current offset
            for (size_t i = 0; i < key.size(); ++i) {
                if (key[i] == lock[(offset + i) % lock.size()]) {
                    current_match++;
                }
            }
            bestMatch = std::max(bestMatch, current_match);

            // Exit early if a perfect match is found
            if (bestMatch == static_cast<int>(key.size())) {
                break;
            }
        }
        return static_cast<double>(bestMatch) / static_cast<double>(key.size());  // Return the ratio of matched bits in best match to length of key
    }

    // IVY:
    // New mutate which takes in a uint64_t genome instead of vector of ints
    static uint64_t mutate(const uint64_t genome, double mutationRate, double indelRate, double sizeRate, int minLen, int maxLen) {
        uint64_t newGenome = genome;
        int gs = calculateDataLength(genome);

        // Flip bits mutation
        int numMutations = Random::getBinomial(gs, mutationRate);
        std::vector<int> temp(gs); // used to keep track of which locations have not been mutated
        std::iota(temp.begin(), temp.end(), 0);

        for (int i = 0; i < numMutations; i++) {
            int pickIndex = Random::getIndex(temp.size());
            newGenome ^= (1 << pickIndex);
            temp[pickIndex] = temp.back();
            temp.pop_back();
        }

        // copy chunk mutation
        if (Random::P(indelRate)) {
            // You only want to take a small chunk of the genome and ensure
            // no copying off the end shenanigans
            // IVY NOTE: This should work (I tested it), the only thing is whether copySize is able
            // to go off the end still... that MIGHT break things
            // IVY NOTE: Couldn't think of any edge cases (unless genome gets too small?), but this
            // should probably get checked. The logic after copysize should be fine
            int copySize = std::max(1, Random::getInt(gs * 0.25, gs * 0.5)) + 1;
            int copyStart = Random::getInt(copySize - 1, gs - 1);
            int destStart = Random::getInt(copySize - 1, gs - 1);
            int diff = copyStart - destStart;

            // Create mask for what elements to take
            uint64_t copyMask = generateMask(copyStart, copySize);
            uint64_t destMask = ~generateMask(destStart, copySize);

            uint64_t extractedChunk = newGenome & copyMask; // Copy chunk into new variable
            newGenome &= destMask;                          // Takes out dest chunk (set to all 0)

            // This means that copy is to the right of dest, so left shift
            if (diff < 0) {
                extractedChunk <<= (diff * -1);
            } else {        // Copy is to the left of the dest, so right shift over to right place
                extractedChunk >>= diff;   
            }

            newGenome |= extractedChunk;                    // Shift chunk over to right place and copy
        }

        // change size mutation
        // IVY NOTE: In the orig. implementation, you only add and take from "end" (I'm considering that to be the right side)
        // This is just the easiest way to do it
        if (sizeRate > 0) {
            if (Random::P(sizeRate)) {
                int loc = Random::getIndex(gs);
                if (Random::P(.5)) {  // Add one to the end
                    if (gs < maxLen) { // unless it's too long
                        gs = (gs << 1) | Random::getInt(0, 1);
                    }
                }
                else {
                    if (gs > minLen) {  // Remove one from the end if not too short
                        gs >>= 1;
                    }
                }
            }
        }
        return newGenome; 
    }

    static std::vector<int> mutate(const std::vector<int>& genome, double mutationRate, double indelRate, int minLen) {
        std::vector<int> newGenome = genome;

        int numMutations = Random::getBinomial(genome.size(), mutationRate);
        std::vector<int> temp(newGenome.size()); // used to keep track of which locations have not been mutated
        std::iota(temp.begin(), temp.end(), 0);

        for (int i = 0; i < numMutations; i++) {
            int pickIndex = Random::getIndex(temp.size());
            newGenome[temp[pickIndex]] = 1 - newGenome[temp[pickIndex]]; // Flip the value
            temp[pickIndex] = temp.back();
            temp.pop_back();
        }

        if (Random::P(indelRate)) {
            int genomeSize = newGenome.size();
            int copySize = std::max(1, Random::getInt(genomeSize * 0.25, genomeSize * 0.5)) + 1;
            int copyStart = Random::getIndex(genomeSize - copySize);
            int dest = Random::getIndex(genomeSize - copySize);

            // Ensure the destination range is within the vector bounds
            if (dest + copySize > genomeSize) {
                dest = genomeSize - copySize; // Adjust destination to fit within bounds
            }

            if ((copyStart < dest && dest < copyStart + copySize) || (dest < copyStart && copyStart < dest + copySize)) {
                // If there's overlap, use std::copy_backward to prevent overwriting source data
                std::copy_backward(newGenome.begin() + copyStart, newGenome.begin() + copyStart + copySize, newGenome.begin() + dest + copySize);
            }
            else {
                // No overlap or safe to use std::copy
                std::copy(newGenome.begin() + copyStart, newGenome.begin() + copyStart + copySize, newGenome.begin() + dest);
            }
        }

        return newGenome;
    }

    static int countMismatchesWithOffsetOLD(const std::vector<int>& longer, const std::vector<int>& shorter, int offset) {
        int mismatches = 0;
        for (int i = 0; i < longer.size(); ++i) {
            // Calculate the corresponding index in the shorter vector, considering the offset
            int shorterIndex = i - offset;

            // If the index is outside the shorter vector, count as a mismatch
            if (shorterIndex < 0 || shorterIndex >= shorter.size()) {
                mismatches++;
            }
            else if (longer[i] != shorter[shorterIndex]) { // Else, compare the values
                mismatches++;
            }
        }
        // Count mismatches for any "overhang" beyond the end of the longer vector
        if (offset + shorter.size() > longer.size()) {
            mismatches += offset + shorter.size() - longer.size();
        }
        return mismatches;
    }

    static int countMismatchesWithOffset(const std::vector<int>& longer, const std::vector<int>& shorter, int offset) {
        int mismatches = 0;
        int longerIndex = offset;
        for (int i = 0; i < shorter.size(); ++i) {
            // Calculate the corresponding index in the shorter vector, considering the offset

            // If the index is outside the shorter vector, count as a mismatch
            if (longer[longerIndex++] != shorter[i]) { // Else, compare the values
                mismatches++;
            }
        }
        return mismatches;
    }

    // IVY: Rewrote to deal with uint64
    static int findMinMismatches(uint64_t A, uint64_t B) {
        int As = calculateDataLength(A);
        int Bs = calculateDataLength(B);

        uint64_t longer  = As >= Bs ? A : B;
        uint64_t shorter = As <= Bs ? A : B;

        // IVY NOTE: Not super sure, what you mean exactly by overhangs, so this might be wrong
        // Basic idea is that you generate masks of right size then shift over and count how many
        // mismatches
        int endOffset = calculateDataLength(longer) - 1;
        int maskSize = calculateDataLength(shorter);
        int minMismatches = INT_MAX;

        // Remove the leading 1 before doing any comparisons
        // NOTE: This is important to do this after you have already calculated mask length
        longer ^= (1 << calculateDataLength(longer));
        shorter ^= (1 << calculateDataLength(shorter));

        for (int offset = maskSize - 1; offset <= endOffset; ++offset) {
            // Generate chunk from longer and line up with shorter
            uint64_t chunk = longer & generateMask(offset, maskSize);
            chunk >>= (offset + 1 - maskSize);
            minMismatches = std::min(minMismatches, __builtin_popcountll(chunk ^ shorter));
        }
        return minMismatches;
    }

    static int findMinMismatches(const std::vector<int>& A, const std::vector<int>& B) {
        const std::vector<int>& longer = A.size() >= B.size() ? A : B;
        const std::vector<int>& shorter = A.size() < B.size() ? A : B;

        int minMismatches = INT_MAX;
        // Iterate over all possible offsets, including full "overhang"
        //for (int offset = -1 * (shorter.size() - 1); offset <= static_cast<int>(longer.size()); ++offset) {
        // Iterate over all possible offsets, no "overhang"
        for (int offset = 0; offset <= static_cast<int>(longer.size()) - static_cast<int>(shorter.size()); offset++) {
            int mismatches = countMismatchesWithOffset(longer, shorter, offset);
            minMismatches = std::min(minMismatches, mismatches);
        }

        return minMismatches;
    }

    // Function to return the number of bits that differ between two unsigned long long numbers
    // IVY: Made more efficient, you can just use built in function to count 1s set
    static uint64_t diffBits(uint64_t a, uint64_t b) {
        auto diff = a ^ b;
        return __builtin_popcountll(diff);

        // IVY: This was your old impl, if you think mine is correct than feel free to remove it
        // unsigned int count = 0;
        // while (diff) {
        //     count += diff & 1;
        //     diff >>= 1; // Shift right by 1 bit
        // }
        // return count;
    }







    // Helper function to calculate the data length(number of bits after the marker bit)
    static int calculateDataLength(uint64_t genome) {
        if (genome == 1) return 0; // Case for where only leading 1 is set
        return 64 - __builtin_clzll(genome) - 1; // -1 to adjust for the leading 1 which is set
    }

    // same function, use this one if there is a compiler problem
    static int calculateDataLength_alt(uint64_t encoded) {
        int length = 0;
        while (encoded >>= 1) {
            length++;
        }
        return length; // Excludes the marker bit
    }

    // Function to find the best match that takes uint64_t
    static std::pair<int, uint64_t> findBestMatch(const uint64_t genome, const uint64_t target) {

        int genomeDataLength = calculateDataLength(genome);
        int targetDataLength = calculateDataLength(target);

        int bestMatchCount = 0;
        uint64_t bestMatchPositions = 0;

        // Slide only within the data part of the genome
        for (int shift = 0; shift <= genomeDataLength - targetDataLength; ++shift) {
            int matchCount = 0;
            uint64_t matchPositions = 0;

            for (int i = 0; i < targetDataLength; ++i) {
                bool genomeBit = (genome >> (genomeDataLength - 1 - shift - i)) & 1;
                bool targetBit = (target >> (targetDataLength - 1 - i)) & 1;

                if (genomeBit == targetBit) {
                    matchCount++;
                    matchPositions |= 1ULL << (genomeDataLength - 1 - shift - i);
                }
            }

            if (matchCount > bestMatchCount) {
                bestMatchCount = matchCount;
                bestMatchPositions = matchPositions;
            }
        }

        return { bestMatchCount, bestMatchPositions };
    }


    // version that takes vector<ints>
    static std::pair<int, uint64_t> findBestMatch(const std::vector<int>& genomeSeq, const std::vector<int>& targetSeq) {
        uint64_t genome = encodeGenome(genomeSeq);
        uint64_t target = encodeGenome(targetSeq);
        return(findBestMatch(genome, target));
    }





    // IVY TODO: Still need to convert this to take in a vector of uint64_t
    static std::vector<double> getFit(const std::vector<std::vector<int>>& pop, const double SHARE, const double FWR, const std::vector<int> FWRsizes, const double TARGET, const std::vector<std::vector<int>> & targets) {

        int popSize = pop.size();

        // count of each genotype
        std::unordered_map<uint64_t, int> uniqueCounts; // count of each genotype
        std::unordered_map<uint64_t, std::vector<int>> reverseLookup; // links each encoded genome to the bit string
        std::vector<uint64_t> encodedPop(pop.size()); // for each member of the pop, store it's encoding


        for (int id = 0; id < popSize; id++) {
            uint64_t encoded = encodeGenome(pop[id]);
            reverseLookup[encoded] = pop[id];
            encodedPop[id] = encoded;
            uniqueCounts[encoded]++;
        }

        std::vector<double> scores(popSize, 0);
        //int complexity_target_max = std::accumulate(targets.begin(), targets.end(), 0,
        //    [](int sum, const std::vector<int>& innerVec) {
        //        return sum + innerVec.size();
        //    });

        std::vector<uint64_t> encodedTargets;
        for (const std::vector<int>& target : targets) {
            encodedTargets.push_back(encodeGenome(target));
        }

        if (TARGET != 0.0) {
            // target coverage just tells us how well each target is being scored - it has no effect on score or modes
            std::vector<double> targetCoverage(targets.size(), 0.0);

            std::unordered_map<uint64_t, double> TARGET_scores;
            std::unordered_map<uint64_t, double> TARGET_complexities_active; // this will track which sites are active at least once

            for (auto iter = uniqueCounts.begin(); iter != uniqueCounts.end(); ++iter) {
                double targetTotal = 0.0;
                int t = 0;
                uint64_t activePositions = 0;
                for (const uint64_t& target : encodedTargets) {
                    auto [bestMatchCount, bestMatchPositions] = findBestMatch(iter->first, target);
                    activePositions |= bestMatchPositions;
                    double targetMatchRate = static_cast<double>(bestMatchCount) / static_cast<double>(targets[t].size());
                    targetTotal += targetMatchRate;
                    targetCoverage[t++] += targetMatchRate * iter->second;
                    if (0) {
                        std::cout << "  " << g2s(reverseLookup[iter->first]) << std::endl << "  " << g2s(decodeGenome(target)) << std::endl
                            << "    " << bestMatchCount << " >> " << targetMatchRate << " || " << targetTotal << std::endl;
                        std::cout << "  matchMask: " << u642s(bestMatchPositions) << std::endl;
                    }
                }
                TARGET_complexities_active[iter->first] = __builtin_popcountll(activePositions);
                TARGET_scores[iter->first] = TARGET * targetTotal / static_cast<double>(targets.size());
                if (0) {
                    std::cout << "    " << TARGET_scores[iter->first] <<
                        "\n  complexity: " << u642s(activePositions) << "  " <<
                        TARGET_complexities_active[iter->first] << std::endl;
                }
            }

            ///////////////////////////////////////////////////////
            // save target coverage
            std::string tempStr = "";
            std::string tempHeaderStr = "";
            for (int t = 0; t < targets.size(); t++) {
                tempHeaderStr += "t" + std::to_string(t) + ",";
                tempStr += std::to_string(targetCoverage[t] / static_cast<double>(popSize)) + ",";
            }
            tempStr.pop_back();
            tempHeaderStr.pop_back();
            FileManager::openAndWriteToFile("targetCoverage.csv",tempStr, tempHeaderStr);
            FileManager::closeFile("targetCoverage.csv");
            // end save target coverage
            ////////////////////////////////////////////////////////

            //std::cout << "TARGET scores: ";
            for (int i = 0; i < popSize; i++) {
                scores[i] += TARGET_scores[encodedPop[i]]; // the more mismathces the worse the score
                //std::cout << TARGET_scores[encodedPop[i]] << " ";
            }
            //std::cout << std::endl;
        }
        
        // Fit When Rare

        if (FWR != 0.0 || SHARE != 0.0) {
            // FWR is optimized to only consider one instance of each uniqe genotype. It thus needs to know
            // the total number of types and the number of each type and it must be able to assigne scores
            // after evaluation to individuals correctly based on their actual type.

            // for each genotype, the diff to each other genotype
            std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>> uniqueDiffs;

            for (auto it1 = uniqueCounts.begin(); it1 != uniqueCounts.end(); ++it1) {
                for (auto it2 = it1; it2 != uniqueCounts.end(); ++it2) {
                    auto diffs = diffBits(it1->first, it2->first);
                    uniqueDiffs[it1->first][it2->first] = diffs;
                    uniqueDiffs[it2->first][it1->first] = diffs;
                }
            }

            // for each genotype, it's FWR_score
            std::unordered_map<uint64_t, double> FWR_scores;

            for (auto it1 = uniqueDiffs.begin(); it1 != uniqueDiffs.end(); ++it1) {
                FWR_scores[it1->first] = 0;
                for (auto it2 = uniqueDiffs[it1->first].begin(); it2 != uniqueDiffs[it1->first].end(); ++it2) {
                    double normalizedFrequency = static_cast<double>(uniqueCounts[it2->first]) / static_cast<double>(popSize); // how common is the genome we are comparing to?
                    FWR_scores[it1->first] += normalizedFrequency
                        * (std::max(std::min(static_cast<double>(FWRsizes[1]), uniqueDiffs[it1->first][it2->first]), static_cast<double>(FWRsizes[0]))
                        -
                            static_cast<double>(FWRsizes[0]));
                    if (0) {
                        std::cout << "  " << u642s(it1->first) << std::endl;
                        std::cout << "  " << u642s(it2->first) << std::endl;
                        std::cout << "     " << static_cast<double>(uniqueCounts[it2->first]) << " " << normalizedFrequency << " : " << uniqueDiffs[it1->first][it2->first] << "  " <<
                            (normalizedFrequency
                            * (std::max(std::min(static_cast<double>(FWRsizes[1]), uniqueDiffs[it1->first][it2->first]), static_cast<double>(FWRsizes[0]))
                                - static_cast<double>(FWRsizes[0]))) << std::endl;
                    }
                }
                FWR_scores[it1->first] /= FWRsizes[1]- FWRsizes[0];
                //std::cout << "  <<<  " << FWR_scores[it1->first] << std::endl;
            }

            //std::cout << "FWR scores: ";
            for (int i = 0; i < popSize; i++) {
                //std::cout << FWR_scores[encodedPop[i]] << " ";
                scores[i] += FWR_scores[encodedPop[i]] * FWR;
                scores[i] += SHARE / ((1.0 - FWR_scores[encodedPop[i]]) * popSize); // use FWR score to figure out the number of genomes like this one and divvy up SHARE
                //std::cout << FWR_scores[encodedPop[i]] << " : " << SHARE << " /  1 - " << FWR_scores[encodedPop[i]]  << " * " << popSize << " = " << SHARE / ((1.0 - FWR_scores[encodedPop[i]]) * popSize) << " || ";
            }
            //std::cout << std::endl;
            //std::cout << std::endl;
        }

        return scores;
    }

    static std::vector<int> roulette_wheel_selection(std::vector<double>& scores, std::vector<double> cumulative_scores, int numPicks = -1) {

        if (numPicks <= 0) {
            numPicks = scores.size();
        }

        std::vector<int> selection_indices;
        if (scores.empty()) {
            return selection_indices; // Return empty if no scores
        }

        // Create a cumulative sum list
        if (cumulative_scores.empty()) {
            cumulative_scores.resize(scores.size());
            cumulative_scores[0] = std::max(0.0, scores[0]);
            for (size_t i = 1; i < scores.size(); ++i) {
                if (scores[i] < 0) scores[i] = 0;
                cumulative_scores[i] = cumulative_scores[i - 1] + scores[i];
            }
        }
        double total_score = cumulative_scores.back();
        selection_indices.reserve(numPicks);
        if (total_score == 0) {
            //std::cout << "ALL DEAD!" << std::endl;
            for (size_t i = 0; i < numPicks; ++i) {
                selection_indices.push_back(Random::getIndex(scores.size()));
            }
        }
        else {
            for (size_t i = 0; i < numPicks; ++i) {
                double pick = Random::getDouble(0, total_score);
                auto it = std::upper_bound(cumulative_scores.begin(), cumulative_scores.end(), pick);
                // std::cout << "choosing: " << std::distance(cumulative_scores.begin(), it) << '\n' << std::flush;
                selection_indices.push_back(std::distance(cumulative_scores.begin(), it));
            }
        }

        return selection_indices;
    }


    // a MODES_TRACKER will manage a compete instance of a tracking system
    // this will also need to keep a current phylo tracker so a function to update this per brith will be needed
    // files saved by this MODES_TRACKER will include the "name", that is, the filterDepth for this tracker
    // I will assume that change delay and filter depth are the same

    // should the trackers update eachothers fridges?

    class MODES_TRACKER {
    public:
        int filterDepth;
        int changeDelay;
        bool saveGenomes;
        int popSize;
        int trackerStart;
        std::string name;
        std::vector<std::vector<int>> popHistory;//(popSize);
        std::vector<int> phylo_current;// (phylo_m_current);
        std::vector<int> phylo_next;// (phylo_m_current);
        std::vector<int> phylo_history;// (phylo_m_current);
        std::vector<double> score_history;//(popSize);

        std::unordered_set<uint64_t> fridge; // longterm storage for novelty metric
        std::queue<std::unordered_set<uint64_t>> fridgeQueue; // stores genomes that have been collected but are not yet in the fridge
    
        MODES_TRACKER(int _filterDepth, bool _saveGenomes, int _popSize, int _trackerStart) :
            filterDepth(_filterDepth),
            changeDelay(_filterDepth),
            saveGenomes(_saveGenomes),
            popSize(_popSize),
            trackerStart(_trackerStart),
            name("filterDepth_" + std::to_string(filterDepth) + "_" + std::to_string(trackerStart))
        {
            popHistory.resize(popSize);
            phylo_current.resize(popSize); // used to track phylogony
            phylo_next.resize(popSize); // used to track phylogony
            std::iota(phylo_current.begin(), phylo_current.end(), 0); // Fill with 0, 1, 2, ...
            phylo_history = phylo_current;
            score_history.resize(popSize, 0);

            if (changeDelay % filterDepth != 0) {
                std::cout << "in tracker " << name << " ::   changeDelay(" << changeDelay << ") % filterDepth(" << filterDepth << ") != 0 - please correct and restart.exitting." << std::endl;
                exit(0);
            }

        }


        void updateMODES(const std::vector<std::vector<int>> & pop, int gen, const std::vector<std::vector<int>> & mPop_locks) {
            if ((gen - trackerStart) % filterDepth == 0) {

                //std::cout << "  Tracker " << name << " :: storing genomes, generation: " << gen << std::endl;

                // first step: identify persistant genomes and save them (if saveGenomes)
                // also, build the ecologyBank (a store of the number of each type of genome in the population, including unfiltered genotypes)
                // and if a genome is persistant, add to the fridge queue
                fridgeQueue.push(std::unordered_set<uint64_t>());
                std::unordered_map<uint64_t, int> ecologyBank; // this will store for every save generation the number of each type in the population. i.e., the "eco counts"
                int persistantCount = 0; // number of genomes that get though the persistance filter

                std::vector<uint64_t> persistentGenomes; // this is the list of genomes that were persitant, used to generate complexities - CAN INCLUDE COPIES of the same genome

                std::string GENOME_DATA = "";
                for (size_t id = 0; id < popSize; ++id) {
                    uint64_t encodedGenome = encodeGenome(popHistory[id]);
                    ecologyBank[encodedGenome]++; // put everything in the ecology bank, we will use fridgeQueue to determine which ones get though the persistance filter later

                    if (std::find(phylo_current.begin(), phylo_current.end(), id) != phylo_current.end()) {
                        if (saveGenomes && gen > 0) {
                            GENOME_DATA += std::to_string(gen - filterDepth) + ",";
                            GENOME_DATA += std::to_string(id) + ",";
                            GENOME_DATA += std::to_string(phylo_history[id]) + ",";
                            GENOME_DATA += std::to_string(score_history[id]) + ",\"";
                            for (int bit : popHistory[id]) {
                                GENOME_DATA += std::to_string(bit);
                            }
                            GENOME_DATA += "\"\n";
                        }
                        // this genome survived... it is important. It will need to go in the fridge!
                        fridgeQueue.back().insert(encodedGenome);
                        persistantCount++;
                        persistentGenomes.push_back(encodedGenome);
                    }
                }

                if (saveGenomes) {
                    //FileManager::openAndWriteToFile("Genomes_" + std::to_string(gen) + "__" + name + ".csv", GENOME_DATA, "generation,id,ancestor,score,genome");
                    //FileManager::closeFile("Genomes_" + std::to_string(gen) + "__" + name + ".csv");
                    FileManager::openAndWriteToFile("Genomes_" + std::to_string(0) + "__" + name + ".csv", GENOME_DATA, "generation,id,ancestor,score,genome");
                    FileManager::closeFile("Genomes_" + std::to_string(0) + "__" + name + ".csv");
                }


                if (fridgeQueue.size() > (changeDelay / filterDepth)) { // we have enough data to start collecting MODES data!

                    double novelCount = 0;
                    double changeCount = 0;

                    fridge.insert(fridgeQueue.front().begin(), fridgeQueue.front().end());

                    // CHANGE and NOVELTY
                    for (const auto& genome : fridgeQueue.back()) {
                        // Check if genome is not in fridge
                        if (fridge.find(genome) == fridge.end()) { // this genome is not in the fridge, so its both novel and a change
                            novelCount++;
                            changeCount++;
                        }
                        else { // maybe change
                            if (fridgeQueue.front().find(genome) == fridgeQueue.front().end()) {
                                changeCount++; // This was seen before (it's in the fridge), but was not present at gen - changeDelay, so it counts as a change
                            }
                        }
                    }


                    // we need to get a count of the number of genotypes that match genotypes that got though the filter for ecology
                    // maybe move up
                    double totalSurvivingTypeCount = 0.0;// count of all genomes who had atleast one of their type persist
                    for (size_t id = 0; id < popSize; ++id) {
                        uint64_t encodedGenome = encodeGenome(popHistory[id]); // for each genome in the histoy
                        if (std::find(fridgeQueue.back().begin(), fridgeQueue.back().end(), encodedGenome) != fridgeQueue.back().end()) {
                            totalSurvivingTypeCount++; // if it matches a genome that past the persistance filter, it counts.
                        }
                    }

                    //ECOLOGY (ENTROPY)
                    double ecology = 0.0;
                    double unfilteredEcology = 0.0;

                    for (const std::pair<uint64_t, int>& genome : ecologyBank) {
                        double unfiltered_p = static_cast<double>(genome.second) / static_cast<double>(popSize);
                        unfilteredEcology += unfiltered_p * log2(unfiltered_p);
                        if (std::find(fridgeQueue.back().begin(), fridgeQueue.back().end(), genome.first) != fridgeQueue.back().end()) {
                            double filtered_p = static_cast<double>(genome.second) / totalSurvivingTypeCount;
                            ecology += filtered_p * log2(filtered_p);
                        }
                    }
                    unfilteredEcology *= -1;
                    ecology *= -1;

                    int filteredUniqueCount = fridgeQueue.back().size(); // number of unique genomes that made it though the persistance filter

                    ////////////////
                    // complexity //
                    ////////////////

                    std::unordered_map<uint64_t, int> uniqueCounts; // count of each genotype
                    std::unordered_map<uint64_t, std::vector<int>> reverseLookup; // links each encoded genome to the bit string
                    std::vector<uint64_t> encodedPop(pop.size()); // for each member of the pop, store it's encoding

                    for (int id = 0; id < popSize; id++) {
                        uint64_t encoded = encodeGenome(pop[id]);
                        reverseLookup[encoded] = pop[id];
                        encodedPop[id] = encoded;
                        uniqueCounts[encoded]++;
                    }

                    std::cout << "  Tracker " << name << " :: Saveing MODES data for (delayed) gen: " << (gen - filterDepth) << "   out of " << persistantCount << " persisten genomes, with " << fridgeQueue.back().size() << " unique types" << std::endl;

                    ///////////////////////
                    // TARGET complexity //
                    ///////////////////////


                    std::vector<uint64_t> encodedTargets;
                    for (const std::vector<int>& target : targets) {
                        encodedTargets.push_back(encodeGenome(target));
                    }

                    uint64_t allActivePositions = 0;
                    double TARGET_complexity = 0;

                    //for (auto iter = uniqueCounts.begin(); iter != uniqueCounts.end(); ++iter) {
                    for (const auto& genome : persistentGenomes) {
                        double targetTotal = 0.0;
                        int t = 0;
                        uint64_t activePositions = 0;
                        for (const uint64_t& target : encodedTargets) {
                            auto [bestMatchCount, bestMatchPositions] = findBestMatch(genome, target);
                            activePositions |= bestMatchPositions;
                            //std::cout << "genome:  " << u642s(genome) << "\ntarget:  " << u642s(target) << "\nmatch:   " << u642s(bestMatchPositions) << "\n         " << u642s(activePositions) << std::endl << std::endl;;
                        }
                        allActivePositions |= activePositions;
                        TARGET_complexity += static_cast<double>(__builtin_popcountll(activePositions))/ static_cast<double>(persistantCount * pop[0].size());
                        //std::cout << __builtin_popcountll(activePositions) << "  " << persistantCount << " " << pop[0].size() << " -> " << TARGET_complexity << std::endl;
                    }                    
                    //std::cout << gen << " : " << name << "  -----------------------  " << TARGET_complexity << std::endl;

                    double all_TARGET_complexity = static_cast<double>(__builtin_popcountll(allActivePositions)) / static_cast<double>(pop[0].size());

                    ////////////////////
                    // FWR complexity //
                    ////////////////////

                    allActivePositions = 0;
                    double FWR_complexity = 0;

                    uint64_t bitmask = (1ULL << pop[0].size()) - 1ULL;

                    for (const auto& genome : persistentGenomes) {
                        uint64_t activePositions = 0;
                        for (const auto& other : uniqueCounts) {
                            if (genome != other.first) {
                                auto [bestMatchCount, bestMatchPositions] = findBestMatch(genome, other.first);
                                activePositions |= ~bestMatchPositions & bitmask;
                                //std::cout << u642s(genome) << "\n" << u642s(other.first) << "\n" << u642s(bitmask) << "\n" << u642s(bestMatchPositions) << "\n" << u642s(activePositions) << std::endl << std::endl;;
                            }
                        }
                        allActivePositions |= activePositions;
                        FWR_complexity += static_cast<double>(__builtin_popcountll(activePositions)) / static_cast<double>(persistantCount * pop[0].size());
                        //std::cout << __builtin_popcountll(activePositions) << "  " << persistantCount << " " << pop[0].size() << " -> " << FWR_complexity << std::endl;
                    }
                    //std::cout << gen << " : " << name << "  -----------------------  " << FWR_complexity << std::endl;

                    double all_FWR_complexity = static_cast<double>(__builtin_popcountll(allActivePositions)) / static_cast<double>(pop[0].size());

                    ////////////////////
                    // SEX complexity //
                    ////////////////////

                    allActivePositions = 0;
                    double SEX_complexity = 0;

                    for (const auto& genome : persistentGenomes) {
                        uint64_t activePositions = 0;
                        for (const auto& sexGenome : mPop_locks) {
                            auto [bestMatchCount, bestMatchPositions] = findBestMatch(genome, encodeGenome(sexGenome));
                            activePositions |= bestMatchPositions;
                            //std::cout << "genome:  " << u642s(genome) << "\nsex gnm: " << u642s(encodeGenome(sexGenome)) << "\nmatch:   " << u642s(bestMatchPositions) << "\n         " << u642s(activePositions) << std::endl << std::endl;;
                        }
                        allActivePositions |= activePositions;
                        SEX_complexity += static_cast<double>(__builtin_popcountll(activePositions)) / static_cast<double>(persistantCount * pop[0].size());
                        //std::cout << __builtin_popcountll(activePositions) << "  " << persistantCount << " " << pop[0].size() << " -> " << SEX_complexity << std::endl;
                    }
                    //std::cout << gen << " : " << name << "  -----------------------  " << SEX_complexity << std::endl;

                    double all_SEX_complexity = static_cast<double>(__builtin_popcountll(allActivePositions)) / static_cast<double>(pop[0].size());

                    /////////////////////
                    // complexity done //
                    /////////////////////


                    std::string MODES_DATA = "";
                    MODES_DATA += std::to_string(gen - filterDepth) + ",";
                    MODES_DATA += std::to_string(changeCount) + ",";
                    MODES_DATA += std::to_string(novelCount) + ",";
                    MODES_DATA += std::to_string(ecology) + ",";

                    MODES_DATA += std::to_string(TARGET_complexity) + ",";
                    MODES_DATA += std::to_string(FWR_complexity) + ",";
                    MODES_DATA += std::to_string(SEX_complexity) + ",";

                    MODES_DATA += std::to_string(all_TARGET_complexity) + ",";
                    MODES_DATA += std::to_string(all_FWR_complexity) + ",";
                    MODES_DATA += std::to_string(all_SEX_complexity) + ",";

                    MODES_DATA += std::to_string(unfilteredEcology) + ",";
                    MODES_DATA += std::to_string(persistantCount) + ",";
                    MODES_DATA += std::to_string(filteredUniqueCount) + ",";
                    MODES_DATA += std::to_string(ecologyBank.size());

                    FileManager::openAndWriteToFile("MODES_data__"+name+".csv", MODES_DATA, "generation,change,novelty,ecology,target_compexity,fwr_complexity,sex_complexity,all_TARGET_complexity,all_FWR_complexity,all_SEX_complexity,unfilteredEcology,filteredCount,filteredUniuqeTypesCount,unfilteredUniqueTypesCount");
                    FileManager::closeFile("MODES_data__" + name + ".csv");

                    fridgeQueue.pop();

                }


                // update history populations and phylogony
                popHistory = pop;
                phylo_history = phylo_current;
                // set each ancestor to self so that the next history will point to the correct ansectors in this current generation
                std::iota(phylo_current.begin(), phylo_current.end(), 0);

            }        
        }

        void updateNextPhylo(int parentID, int offspringID) {
            phylo_next[offspringID] = phylo_current[parentID];
        }

        void updatePhylo() {
            phylo_current = phylo_next;
        }

        void updateScoreHistory(std::vector<double> scores,int gen) {
            // only if on update
            if ((gen - trackerStart) % filterDepth == 0) {
                score_history = scores;
            }
        }
    };



    // parameters for group and brain namespaces
    
    // a local variable used for faster access to the ParameterLink value
        
    class Parasite {
    public:
        std::vector<int> genome;
        double resource;
        int age;

        Parasite(const std::vector<int>& genome, double resource, int age)
            : genome(genome), resource(resource), age(age) {}
    };


    MODESWorld(std::shared_ptr<ParametersTable> PT);
    virtual ~MODESWorld() = default;

    virtual auto evaluate(std::map<std::string, std::shared_ptr<Group>>& /*groups*/, int /*analyze*/, int /*visualize*/, int /*debug*/) -> void override;

    virtual auto requiredGroups() -> std::unordered_map<std::string, std::unordered_set<std::string>> override;
};

