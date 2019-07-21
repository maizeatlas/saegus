__author__ = 'DoubleDanks'
import numpy as np

def test_extract_genotype_matrix(pop, extracted_genotype_matrix, test_result_filename):
    """
    A test for the correctness of the results of the extract_genotype_matrix function.
    Compares the genotypes in the matrix with the genotypes stored in the population.
    Result is written to a json file.
    If the test passes then every entry sums to ploidy * number_of_loci and a .json file
    is written.
    """
    test_results = {}
    pop_ploidy = pop.ploidy()
    ploidy_times_loci = pop.ploidy() * pop.totNumLoci()
    for idx, ind in enumerate(pop.individuals()):
        # using a list comprehension to avoid a nested loop
        result_set = [sum(np.equal(extracted_genotype_matrix[p][idx], ind.genotype(ploidy=p))) for p in range(pop_ploidy)]
        test_results[idx] = result_set
        assert sum(result_set) == ploidy_times_loci
        "Value discovered which is not equal to the expected value."
    result_array = np.array(list(test_results.values()))
    np.savetxt(test_result_filename, result_array, delimiter='\t', fmt='%d')
    return test_results

def test_sampling_of_meta_population(cumulative_metapop_filename, list_of_sampled_filenames, test_result_filename):
    """
    A test to determine that the operator ReplicateMetaPopulation is sampled individuals correctly.
    TestReplicateMetaPopulation is the ReplicateMetaPopulation with the additional step of writing
    .pop files of each sample during run-time.

    We compare the genotypes of the cumulative metapopulation and the corresponding sampled population.
    The expected result is ploidy * num_loci. A file of the results is written for shareability of results.
    :param cumulative_metapop_filename:
    :type cumulative_metapop_filename:
    :param list_of_sampled_filenames:
    :type list_of_sampled_filenames:
    :param test_result_filename:
    :type test_result_filename:
    :return:
    :rtype:
    """
    pass



def test_edited_genotype_to_array_conversion(first_homologous_set, second_homologous_set, complete_geno):
    """
    Test to determine whether the original string, unsplit genotypes can be reassembled from the arrays
    they are split into. Reassemble each genotype and sum the boolean values. If correct then value will
    be equal to number of loci. If not correct value will be less than number of loci.

    :param first_homologous_set:
    :type first_homologous_set:
    :param second_homologous_set:
    :type second_homologous_set:
    :param complete_geno:
    :type complete_geno:
    :return:
    :rtype:
    """
    failing_test_results = {}
    passing_test_results = {}
    for ind in range(105):
        # ind does not change from zero so all genotypes are compared against the first genotype
        failing_test_results[ind] = sum(['/'.join((str(first_homologous_set[0, i]),
                                                   str(second_homologous_set[0, i]))) == complete_geno[ind, i]
                                         for i in range(44445)])
        passing_test_results[ind] = sum(['/'.join((str(first_homologous_set[ind, i]),
                                                   str(second_homologous_set[ind, i]))) == complete_geno[ind, i]
                                         for i in range(44445)])
        test_results_array = np.array([list(failing_test_results.keys()), list(failing_test_results.values()),
                                       list(passing_test_results.values())]).T
    np.savetxt('test_results_genotype_to_array_conversion.txt', test_results_array, fmt='%d', delimiter='\t',
               header='tested by reassembling original genotypes from split arrays and comparing')