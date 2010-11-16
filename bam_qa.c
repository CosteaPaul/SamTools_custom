#include <stdio.h>
#include "sam.h"

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core)
{

  if (core->flag&BAM_FUNMAP)
    return false;
  
  return true;
}

/**
 * Main of app
 */
int main_qa(int argc, char *argv[])
{
  samfile_t *fp;
  FILE *outputFile;
  if (argc == 1) {
    fprintf(stderr, "Usage: qaCompute <in.bam> <output.out>\n Note that the .bam file should be ordered\n");
    return 1;
  }
  //Note that file is supposed to have been ordered beforehand!
  if ((fp = samopen(argv[1], "rb", 0)) == 0) {
    fprintf(stderr, "qaCompute: Fail to open BAM file %s\n", argv[1]);
    return 1;
  }
  if ((outputFile = fopen(argv[2], "wt")) == 0) {
    fprintf(stderr, "qcCompute: Filed to create output file %s\n", argv[2]);
    return 1;
  }
  
    //Initialize bam entity
    bam1_t *b = bam_init1();

    //All var declarations
    int64_t totalGenomeLength = 0;
    int32_t unmappedReads = 0;
    int32_t zeroQualityReads = 0;
    int32_t totalNumberOfReads = 0;
    int32_t totalProperPaires = 0;
    int32_t chrSize = 0;
    int32_t maxCoverage = 30; //Should read from user.
    int *entireChr = NULL;
    //Keep header for further reference
    bam_header_t* head = fp->header;
    
    int32_t currentTid = -1;

    //Create "map" vector for histogram
    int64_t* coverageHist = (int64_t*)malloc((maxCoverage+1)*sizeof(int64_t)); 
    memset( coverageHist, 0, (maxCoverage+1)*sizeof(int64_t));

    //Write file table header
    fprintf(outputFile, "Chromosome\tCoverage\n");

    while (samread(fp, b) >= 0) {
      //Get bam core.
      const bam1_core_t *core = &b->core;
      if (core == NULL) {
	//There is something wrong with the read/file
	printf("Input file is corrupt!");
	//Leak everything and exit!
	return -1;
      }
      //BAM block has been read
      if (!is_mapped(core))
	++unmappedReads;
      else {

	if (core->tid != currentTid) {
	  
	  //Count coverage!
	  if (currentTid != -1) {
	    int32_t covVal = 0;
	    int32_t covSum = 0;
	    int32_t i;

	    //Go through chromosome and count avarage covarage.                                                                                                                                                                                         
	    for (i=0; i<chrSize; ++i){
	      covVal += entireChr[i];
	      covSum += covVal;
	      //Add value to histogram
	      if (covVal > maxCoverage) {
		++coverageHist[maxCoverage];
	      } else {
		++coverageHist[covVal];
	      }

	    }
	    
	    //Printout avarage coverage over this chrom
	    printf("Average cover over %s : %3.2f\n", head->target_name[currentTid], (double)covSum / chrSize);
	    if (currentTid < 24)
	      //Don't print alternates to the file
	      fprintf(outputFile, "%s\t%3.2f\n", head->target_name[currentTid], (double)covSum / chrSize);

	  }

	  //Get length of next section                                                                                       
          chrSize = head->target_len[core->tid];
          totalGenomeLength += chrSize;
          printf("Computing %s... \n",head->target_name[core->tid]);

	  //Done with current section.
	  //Allocate memory
	  entireChr = (int*)realloc(entireChr, (chrSize+5)*sizeof(int));
	  
	  if (entireChr == NULL) {
	    printf("Allocation failed! \n");
	    return -1;
	  }
	  memset(entireChr, 0, (chrSize+5)*sizeof(int));
	  
	  currentTid = core->tid;
	
	}
	
	//If read has quality == 0, we won't count it as mapped
	if (core->qual != 0) {
	  if (core->flag&BAM_FPROPER_PAIR) {
	    //Is part of a proper pair
	    ++totalProperPaires;
	  }
	  //All entries in SAM file are represented on the forward strand! (See specs of SAM format for details)
	  ++entireChr[core->pos];
	  --entireChr[core->pos+core->l_qseq];

	} else {
	  //Count is as unmapped?
	  ++zeroQualityReads;
	}
      }

      ++totalNumberOfReads;
      
    }

    bam_destroy1(b);
    free(entireChr);

    //Print header for next table in output file
    fprintf(outputFile,"\nCov*X\tProcentage\n");

    printf("Total genome lenght %ld \n", totalGenomeLength);
    //Compute procentages of genome cover!.
    int i=0;
    for (i; i<=maxCoverage; ++i) {
      if (i == 0) {
	//Non-covered!
	printf("%3.2f of genome has not been covered\n", (double)(100*coverageHist[i])/totalGenomeLength);
      } else {
	int64_t coverage = 0;
	//All that has been covered i, had been covered i+1, i+2 and so on times. Thus, do this addition
	int x = i;
	for (x; x<=maxCoverage; ++x) coverage += coverageHist[x];
	printf("%3.2f of genome has been covered at least %dX \n", (double)(100*coverage)/totalGenomeLength, i);
	fprintf(outputFile,"%d\t%3.2f\n",i, (double)(100*coverage)/totalGenomeLength);
      }

    }

    fprintf(outputFile,"\nOther\n");

    //Printout procentage of mapped/unmapped reads                                                                                                     
    double procentageOfUnmapped = (100*unmappedReads)/totalNumberOfReads;
    double procentageOfZeroQuality = (100*zeroQualityReads)/totalNumberOfReads;
    fprintf(outputFile,"Total number of reads: %d\n", totalNumberOfReads);
    fprintf(outputFile,"Procentage of unmapped reads: %3.2f\n", procentageOfUnmapped);
    fprintf(outputFile,"Procentage of zero quality mappings: %3.2f\n", procentageOfZeroQuality);
    int32_t nrOfPaires = totalNumberOfReads/2;
    double procOfProperPaires = (double)(100*(double)totalProperPaires/2)/nrOfPaires;
    fprintf(outputFile,"Number of proper paired reads: %d\n", totalProperPaires);
    fprintf(outputFile,"Procentage of proper paires: %3.2f\n", procOfProperPaires);

    printf("Out of %d reads, you have %3.2f unmapped reads\n and %3.2f zero quality mappings", totalNumberOfReads ,procentageOfUnmapped, procentageOfZeroQuality);
    

    free(coverageHist);

  
  samclose(fp);
  fclose(outputFile);
  return 0;
}
