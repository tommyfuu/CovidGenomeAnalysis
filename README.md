# CovidGenomeAnalysis

Covid genome analysis project at the HMC BioMakerspace. Using machine learning to investigate the pattern between location and covid sequences.

### The GISAIDExtract Directory

The working python files are ORFExtractor.py and scaleORFExtract.py, each responsible for different things.

**ORFExtractor.py** is where we can take in an input fasta file (example file see 'GISAIDExtract/EPI_ISL_402124.fasta') and output the rough orf and alignment score with the reference file.
**scaleORFExtract.py** is where we deal with the giant fasta file. You might have downloaded a giant fasta file at a later day than me. For me, I have it as 'GISAIDExtract/sequences_2020-10-29_07-32.fasta', which I stored at the top of the file as a variable **BigFastaAddress** - you can change it as you need. This file contains one function **scaleOrfAlignScore(BigFastaAddress, scoreDict={})** that takes in the big fasta file and an empty dictionary of scores. The keys of the scoreDict will be the orfs available, and the values attached to each key will be all such orfs in the sequences that we looped through and the alignment scores achieved by aligning this orf with the reference rf, and the location and the ascension number this orf is from. An example of this is shown below (with the first 5 sequences in the big fasta file):

```
{'ORF1a': [('IVDC-HB-01', 13230.0, 'Wuhan'), ('IVDC-HB-04', 13230.0, 'Wuhan'), ('IVDC-HB-05', 13230.0, 'Wuhan'), ('IPBCAMS-WH-01', 13227.0, 'Wuhan'), ('WIV04', 13230.0, 'Wuhan')], 'ORF1b': [('IVDC-HB-01', 8106.0, 'Wuhan'), ('IVDC-HB-04', 8106.0, 'Wuhan'), ('IVDC-HB-05', 8104.0, 'Wuhan'), ('IPBCAMS-WH-01', 8106.0, 'Wuhan'), ('WIV04', 8106.0, 'Wuhan')], 'ORFS': [('IVDC-HB-01', 3882.0, 'Wuhan'), ('IVDC-HB-04', 3882.0, 'Wuhan'), ('IVDC-HB-05', 3882.0, 'Wuhan'), ('IPBCAMS-WH-01', 3882.0, 'Wuhan'), ('WIV04', 3882.0, 'Wuhan')], 'ORF3a': [('IVDC-HB-01', 873.0, 'Wuhan'), ('IVDC-HB-04', 873.0, 'Wuhan'), ('IVDC-HB-05', 873.0, 'Wuhan'), ('IPBCAMS-WH-01', 873.0, 'Wuhan'), ('WIV04', 873.0, 'Wuhan')], 'ORFN': [('IVDC-HB-01', 1302.0, 'Wuhan'), ('IVDC-HB-04', 1301.0, 'Wuhan'), ('IVDC-HB-05', 1302.0, 'Wuhan'), ('IPBCAMS-WH-01', 1302.0, 'Wuhan'), ('WIV04', 1302.0, 'Wuhan')]}
```

From this dictionary, we extract information for each orf - for each orf, we make a pandas df and a corresponding csv file for that csv, containing all the ascensionNumbers, alignment scores, and the locations. For example, check out **ORF1aAlignmentScore.csv**. The current files are just from the first 5 sequences in the big fasta file. To process all the files and generate comprehensive csv files, I am currently running our functions with the help of google cloud computing services. They will be uploaded once this is done.

### Next steps

1. When we download a new giant files, we decide what are new in the new giant fasta, and then update the scoreDict and csv files accordingly
2. Possibly automate step 1
3. Normalize alignment scores
4. Building clustering model (pair programming!)
5. Possibly incorporating muscle into our project if popssible
6. Look for better orf prediction methods and incorporate that into our model

### Acknowledgement

Developers:\
Tom Fu [(@tommyfuu)](https://github.com/tommyfuu)\
Lucy Paddock [(@lucinda-paddock)](https://github.com/lucinda-paddock)\
Liam Chalk [(@liamchalk00)](https://github.com/liamchalk00)\
April Zhao

Affiliation: Harvey Mudd College BioMakerspace (Polymerspace). \
Faculty advisor: Dr. Dan Stoebel.

This project is part of the intellectual property of the [Harvey Mudd College BioMakerspace](https://biomakerspace.com/). Please notify the space at **tfu@g.hmc.edu** for consent before using any of the materials in this repository. Please cite us if you use our repository for academic sharing purposes.

The funding support of the Harvey Mudd College BioMakerspace comes from Harvey Mudd College Office of Community Engagement and the college's Shanahan Fund.
