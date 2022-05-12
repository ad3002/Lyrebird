#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 25.12.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
from PyExp import AbstractModel

class OrganismModel(AbstractModel):
	''' Container for information about taxon.

	Dumpable attributes:

	- organism_taxon
	- organism_common_name
	- organism_acronym
	- organism_description
	- organism_wgs_projects
	- organism_genome_assemblies

	'''

	dumpable_attributes = [
            "organism_taxon",
            "organism_common_name",
            "organism_acronym",
            "organism_description",
            "organism_wgs_projects",
            "organism_genome_assemblies",
            ]

