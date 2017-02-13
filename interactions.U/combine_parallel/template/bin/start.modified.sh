#***************************************************************************
#                             catRAPID 
#                             -------------------
#    begin                : Nov 2010
#    copyright            : (C) 2010 by G.G.Tartaglia
#    author               : : Tartaglia
#    date                 : :2010-11-25
#    id                   : caRAPID protein-RNA interactions
#    email                : gian.tartaglia@crg.es
#**************************************************************************/
#/***************************************************************************
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU General Public License as published by  *
# *   the Free Software Foundation; either version 2 of the License, or     *
# *   (at your option) any later version.                                   *
# *                                                                         *
# ***************************************************************************/
# please read and install tester.cpp

#!/bin/bash
# takes the sequences (label - sequence)
#echo  $1 $2 > ./database/protein.txt

#cp  ./database/rna.50.3000.txt  ./database/rna.txt
  
awk '{print $1,"'$1'"}' ./protein.list.txt | sed 's/>//g' > ./database/interactions.protein.50.750.txt 

# computes the score
./bin/tester 1 1 50  >    tmp/log.tmp

# calculates the score
bash ./bin/processer.sh positives.txt > tmp/score.txt

# The last column is the size of the RNA
paste ./database/interactions.protein.50.750.txt tmp/score.txt | awk '{print $0}' | sed 's/ENSU/ENST/g'

