
#include "fenrichcpp.h"

using namespace std;

void fenrich_cpp::readGenotypes(string fvcf){
    // Opening files
    bcf_srs_t * sr = bcf_sr_init();
    if(!(bcf_sr_add_reader (sr, fvcf.c_str()))) {
		switch (sr->errnum) {
		case not_bgzf: cout << "File not compressed with bgzip!" << endl; break;
		case idx_load_failed: cout << "Impossible to load index file!"<< endl; break;
		case file_type_error: cout << "File format not detected by htslib!" << endl; break;
		default : cout << "Unknown error!" << endl;
		}
	}
    int n_samples = bcf_hdr_nsamples(sr->readers[0].header);

    unsigned int linecount=0;
    // Read genotype data 
    int ngt, ngt_arr = 0, nds, nds_arr = 0, * gt_arr = NULL, nsl, nsl_arr = 0, * sl_arr = NULL;
    float * ds_arr = NULL;
    bcf1_t * line;
    while(bcf_sr_next_line (sr)){
        linecount ++;
        if (linecount % 100000 == 0) cout << "Read " << to_string(linecount) << " lines" << endl;
        line = bcf_sr_get_line(sr,0);
        if(line->n_allele == 2){
            ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);
            if(ngt == 2*n_samples){
                string sid = string(line->d.id);
                string chr = string(bcf_hdr_id2name(sr->readers[0].header, line->rid));
                int pos = line->pos +1;

                genotype_id.push_back(sid); // Read variant ID
                genotype_chr.push_back(chr); // Read variant chr
                string genotype_ref = string(line->d.allele[0]); // Read reference allele.
                genotype_start.push_back(pos);
                nsl = bcf_get_info_int32(sr->readers[0].header, line, "END", &sl_arr, &nsl_arr);
                if (nsl >= 0 && nsl_arr == 1) genotype_end.push_back(sl_arr[0]);
                else genotype_end.push_back(genotype_start.back() + genotype_ref.size() - 1);
                genotype_val.push_back(vector < float > (n_samples, 0.0));

                for(int i = 0; i < n_samples ; i ++) {
                    if (gt_arr[2*i+0] == bcf_gt_missing || gt_arr[2*i+1] == bcf_gt_missing) bcf_float_set_missing(genotype_val.back()[i]);
                    else genotype_val.back()[i] = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]);
				}
		    } 
        }
    }
    genotype_count = linecount;
    cout << "Read " << to_string(genotype_count) << endl;
}