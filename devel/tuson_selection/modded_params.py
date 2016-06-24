def randomly_convert_fixed_sites(pop, fixed_sites)
    
    Randomly converts fixed sites in pop to nonfixed_sites by changing the allele state at that site.
    
    alleles = [0, 1, 2, 3]
    random.shuffle(fixed_sites)
    for site in fixed_sites
        random_id = random.choice(pop.indInfo('ind_id'))
        random_individual = pop.indByID(random_id)
        current_allele_state = random_individual.allele(site)
        possible_replacements = [allele for allele in alleles if
                                 allele != current_allele_state]
        replacement_allele = random.choice(possible_replacements)
        random_individual.setAllele(replacement_allele, site)