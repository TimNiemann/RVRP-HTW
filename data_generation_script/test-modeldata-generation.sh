#!/bin/bash

TESTDATA="../test-data"
TESTFOLDERS=('paper-evaluation-small-tw')
CITIES=('bayern' 'berlin' 'brandenburg' 'essen' 'hamburg' 'herne' 'hessen' 'muenchen')

# Generating Instances
cd ../instance_generator
make
for i in {0..4}
do
	for folder in "${TESTFOLDERS[@]}"
	do
		echo "Generating ${folder}-${i}"
		mkdir -p "$folder-${i}"
		for city in "${CITIES[@]}"
		do
			#echo "${city}"
			echo "./bin/instanceGenerator ${TESTDATA}/paper-evaluation/${folder}-${i} ${city} -o ${folder}-${i}/${city}"
			./bin/instanceGenerator "${TESTDATA}/paper-evaluation/${folder}-${i}" "${city}" -o "${folder}-${i}/${city}"
		done
	done
done

# Preprocessing Instances
cd ../preprocessing
make
for i in {0..4}
do
	for folder in "${TESTFOLDERS[@]}"
	do
		echo "Preprocessing ${folder}-${i}"
		for city in "${CITIES[@]}"
		do
			#echo "${city}"
			mkdir -p ${TESTDATA}/model-data-paper/${folder}-${i}/${city}
			for customers in 20 40 60 80 100 200 300 400 500
			do
				for days in 5 10 15
				do
					for windows in 2 4 6 8
					do
						echo "./bin/preprocessor ../instance_generator/${folder}-${i}/${city} -r ${customers} -d ${days} -w 0.${windows} -o ${TESTDATA}/model-data-paper/${folder}-${i}/${city}/${city}${i}_r${customers}_d${days}_w0.${windows}.dat"
						./bin/preprocessor "../instance_generator/${folder}-${i}/${city}" -r "${customers}" -d "${days}" -w "0.${windows}" -o "${TESTDATA}/model-data-paper/${folder}-${i}/${city}/${city}${i}_r${customers}_d${days}_w0.${windows}.dat"
					done
				done
			done			
		done
	done
done

# Cleaning Up
cd ../instance_generator
for i in {0..4}
do
	for folder in "${TESTFOLDERS[@]}"
	do
		echo "Cleaning ${folder}-${i}"
		rm -r "$folder-${i}"
	done
done

echo "Done."
