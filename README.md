1)  Creating GlyGen Datasets (dataset-maker directory)
    The following process downloaded data (https://data.glygen.org/ln2downloads/) based
    on dataset master list which is also released with each GlyGen relase. 
    (https://data.glygen.org/ln2data/releases/data/v-1.12.1/misc/dataset-masterlist.json). 

    1) Input and output folders

        These step expects the following symbolic link to input folders that 
        will be holding input data.
            - downloads (https://data.glygen.org/ln2downloads/)
            - compiled  (https://data.glygen.org/ln2data/releases/data/v-1.12.1/compiled)
    
        For output data, create the following folders and make symbolic link to them 
        within the dataset-maker folder.
            - unreviewed
            - reviewed
            - alignments
            - logs


    2) Naming of dataset files
        Species specific datasets are named as $species_$molecule_*.* while those that
        are species agnostic are named as $molecule_*.* where $molecule is one of
        [protein, glycan, proteoform]

    
    3) Making dataset files

        a)  Making protein datasets for all species (*_protein_*.*)
            $ python wrap-make-all-species-mol-specific.py -m protein 

        b)  Making species agnostic glycan datasets (glycan_*.*)
            $ python wrap-make-all-agnostic-mol-specific.py -m glycan
    
        c)  Making species agnostic protein datasets (protein_*.*)
            $ python wrap-make-all-agnostic-mol-specific.py -m protein

        d)  Making proteoform datasets for all species (*_proteoform_*.*)
            $ python wrap-make-all-species-mol-specific.py -m proteoform

        e)  Checking successfull creating of datasets
            $ python wrap-check-datasets-all.py -m protein
            $ python wrap-check-datasets-all.py -m glycan
            $ python wrap-check-datasets-all.py -m protoeoform

        f)  QC checking
            $ python3 qc/check-csv-sanity.py
            $ cat qc/logs/sanity_qc.json | grep file | sort -u



Making data release package

    1)  Making objects for data.glygen.org
    
        a)  Copy dataset files from unreviewed to reviewed
            $ python3 object-maker/copy-dataset-files.py
            $ python3 object-maker/compare-reviewed.py -v $old_ver

        b)  Create BCOs objects in generated/datasets/jsondb/bcodb/
            $ python3 object-maker/make-bcodb.py

        c)  Create BCO extracts in generated/datasets/jsondb/extractdb/
            $ python3 object-maker/make-extractdb.py 

        d)  Create historydb objects 
            $ python3 object-maker/make-historydb-pairs.py 
            $ python3 object-maker/make-historydb-track.py 
    
        e)  Create htmldb objects in generated/datasets/jsondb/htmldb
            $ python3 object-maker/make-htmldb.py
    
        f)  Create initdb object in generated/datasets/jsondb/initdb
            $ python3 object-maker/make-initdb.py 

        g)  Making objects for glygen.org (proteindb/glycandb/ ...)
            $ cd object-maker
            $ sh make-all.sh


    2)  Making data release -- this should be done only from dev server

        a)  Make release dir /data/shared/glygen/releases/data/v-$ver/ and 
            copy files into it
            $ python3 release-maker/release-frontend-data.py -v $ver
	    $ python3 release-maker/release-dsviewer-data.py -v $ver
	
	The script release-dsviewer-data.py will also perform:
		- slinks for motif_ac (for images)
            	- versionize objects in bcodb,extractdb ...
        
	b)  Make ftp release dir /data/shared/glygen/releases/ftp/v-$ver
            $ python3 release-maker/make-ftp-release.py -v $ver


    3) Making triple store release
        
        a)  Make nt files under generated/sparql/glygen/
            
            $ cd triple-maker 
            $ sh wrap-make-glygen-triples.sh 1.12.3
        

        b)  Package triplestore data release
            
            $ cd release-maker
            $ sh make-triplestore-release.sh



Deploying data package
    
    a) The variable $server should be one of [dev, tst, beta, prd]
        $ sh /data/shared/glygen/deployment/deploy-data.sh $server $ver






