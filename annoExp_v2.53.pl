#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);
use Tabix;
use Getopt::Long;
use File::Path;
use Cwd 'abs_path';
use Text::NSP::Measures::2D::Fisher::twotailed;

#build:2016-11-20
#modified:       2016-11-23 Add TABIX to Get REAVEL Database Score;
#modified:v2.0   2016-12-30 Add impact score calculate 
#modified:v2.1   2016-2-11  Add pedigree and interiant model juge;normalized score;mix with gene socre
#                2017-2-16  Add family options for multifamily pedigree file
#modified:v2.2   2016-2-22  Add Gene function,KO,GO;use dbNSFP3.2_gene file
#                2016-2-22  Lower low impact variant ranksore,including disruptive_inframe_deletion,splice_region_variant 
#modified:v2.3   2017-2-27  Add vcf list read model;
#                2017-2-28  Add multifamily model
#modified:v2.4   2017-8-09  fix bug;modified stat;modified dbsnp
#modified:v2.5   2018-1-18  add allele freq,exon rank,fisher_p
#modified:v2.53   2018-1-18  add mulitple gene


my ($vcf,$vcflist,$interp,$sample,$outdir,$help,$pedgree,$otherinfo,$family,@AllSample_vcf,$dbsnp,$minfreq);
my ($genescore);

GetOptions(
    "ped|pedgree:s"=> \$pedgree,
    "list:s"=>\$vcflist,
    "fam|family:s"=>\$family,
    "hit:s"=> \$interp,
    "help"=> \$help,
    "otherinfo"=> \$otherinfo,
    "s|sample:s"=> \$sample,
    "o|outdir:s"=> \$outdir,
    "dbsnp!"=>\$dbsnp,
    "genescore:s"=>\$genescore,
    "minfreq:s"=>\$minfreq,
);


$vcf=$ARGV[0];
$family ||="";
die(&usage) if(!$vcf && !$vcflist || $help);

sub usage{
    print STDERR<<EOF;
    perl $0 <options> [Anno_vcf]
    options:
          -o     outdir 
          -hit   reohit result file
          -list  vcf list
          -genescore  genescore_file 
          -ped    pedgree_file
          -dbsnp  
          -otherinfo print all sample info no matter wheather in pedigree
          -fam|family family name
          -minfreq  kg_eas_af,exac_eas_af,in_house_af,eg:0.01,0.01,0.1 default:0.05,0.05,0.1
          -help
EOF
exit 1;
}


##dbfile

my $bgzip ||="/ifs1/Reohealth/pub/bio-software/samtools/bgzip";
my $tabix ||="/ifs1/Reohealth/pub/bio-software/samtools/tabix";
#my $dbNSFP ||="/ifs1/Reohealth/pub/database/hg19/dbNSFP/dbNSFP.tsv.gz";
#my $dbSNP ||="/ifs1/Reotech/pub/Database/Human/GATK_bundle/2.8/b37/dbsnp_138.b37.vcf.gz";
#my $dbSNP ||="/ifs1/Reotech/pub/Database/Human/GATK_bundle/dbsnp/dbsnp150/common_all_20170710.vcf.gz";
my $dbSNP ||="/ifs1/Reotech/pub/Database/Human/GATK_bundle/dbsnp/dbsnp150/All_20170710.vcf.gz";
$genescore ||="/ifs1/Reotech/pub/biosoft/DNA/DNASeq_bulid/pipeline/PWGS/database/DCM.finalVSreo.list";
my $dieaseinfo ||="$Bin/../../database/Disease_info/Card_Gene_info";
my $Cardinfo ||="$Bin/../../database/Disease_info/Card_missense_ACMG.tsv.gz";
my $reavel ||="/ifs1/Reotech/pub/Database/Human/REVEL/revel_all_chromosomes.tsv.gz";
my $mcap ||="/ifs1/Reotech/pub/Database/Human/MCAP/mcap_v1_0.txt.gz";
my $dbNSFP_gene ||="/ifs1/Reohealth/pub/database/GRCh38/dbNSFP/dbNSFP3.3a/dbNSFP3.3_gene.complete.gz";


my %RankScore;
my %dieaseinfo;
my %Cardinfo;
my %hash;
my %genescore;
my %hNSFP_Gene;
my %standard_score;
my %pedigree;
my %fam_n;
my %group;
my %sample_fam;

%RankScore = (
       'stop_gained'=>10,  
       'stop_lost'=>10,
       'start_lost'=>10,
       'splice_acceptor_variant'=>9,
       'splice_donor_variant'=>9,
       'frameshift_variant'=>9,
       'missense_variant'=>8,
       'disruptive_inframe_deletion'=>4,
       'splice_region_variant'=>4,
       'intron_variant'=>3,
       '3_prime_UTR_variant'=>3,      
       'inframe_deletion'=>3,
       'inframe_insertion'=>3, 
       '5_prime_UTR_premature_start_codon_gain_variant'=>2.5,
       '5_prime_UTR_variant'=>2.5,
       'upstream_gene_variant'=>2.5,
       'downstream_gene_variant'=>2.5,
       'TF_binding_site_variant'=>2,
       'non_coding_exon_variant'=>2,
       'synonymous_variant'=>1.5,
       'intergenic_region'=>1,
       'intragenic_variant'=>1,
       'sequence_feature'=>1,
);

%standard_score = (
    'PVS'=>59.71,
    'PS'=>40.68,
    'PM'=>20.15,
    'PP'=>10.44,
    'BS'=>-40.68,
    'BM'=>-20.15,
    'BP'=>-10.44,
);


$minfreq ||= "0.05,0.05,0.1";
my($kg_eas_af_min,$exac_eas_af_min,$in_house_af_min) = split(/\,/,$minfreq);


#readDieseaInfo($dieaseinfo);
#readCardInfo($Cardinfo);

#globle
my $na_mark=".";
my ($nstopg,$nstopl,$nstartl,$nsplicea,$ndbsnp,$nnovel,$nspliced,$nframe,$inframe_ins,$inframe_del,$inframe_dirup_ins,$inframe_dirup_del,$nintergenic,$nexonice,$nsplicing,$nUTR,$nsplice_region,$nintronic,$nstream,$nnoncoding,$nTFbind,$nintragenic,$ti,$tv,$nvariants,$nmissense,$nsynonymous,$nlof,$nnmd);
my ($name,$path,@out_head);
my $max_impact_score=0;

my @head_cn = qw(染色体 坐标起点 dbsnp_rs Ref Alt 变异类型 纯杂合 基因型别 变异碱基测序深度 基因名称 转录本功能 转录本 
                 危害最大的转录本功能 功能影响分级 功能缺失预测 蛋白截断预测 Mcat有害性预测 Missese多工具综合危害性预测 
                 ExAC东亚人频率 千人东亚人频率 ClinVar临床意义 ClinVar报告 ClinVar级别 HGMD_ID HGMD_目录 HGMD_文献
                 ACMG判断 ACMG证据 疾病 遗传模式 疾病系统 疾病描述);
my @head=qw(SampleName MutType CHROM POS ID  REF ALT Variant_type Impact_score Inheritance_Model Alt_depth 
            Count(case_alt_ref:control_alt_ref) Alt_freq(case:control) fisher_p Gene_symbol_list
            Functional  Triscript High_effect_gene_symbol High_effect High_effect_Triscript Hgv.p Rank_exon Hgv.c Annotation_Impact LOF
            NMD KG_ALL_AF KG_EAS_AF ESP_AF  ExAC_ALL_AF ExAC_EAS_AF In_House_EAP_Freq MCAP_Score  MCAP_pred 
            REVEAL_Score  REAVEL_pred Missense_Multi_pred Affect_splicing Highly_conserved  CGD_Diease  CGD_Inheritance
            CGD_PumMed  OMIM_Disease(gene)  OMIM_Gene ClinVar_ID ClinVar_Clinical_Significance ClinVar_Report  ClinVar_Star  
            ClinVar_Disease HGMD_ID HGMD_Categorization HGMD_Disease  HGMD_PubMed GWAS_OR/Beta  GWAS_trait  
            GWAS_PubMed Repeat_name Functional_domain_name  Category  evidences evidences_info Pathway(KEGG)_id 
            Pathway(KEGG)_full Function_description Disease_description GO_biological_process GO_cellular_component GO_molecular_function);

@out_head=("#$head[0]",@head[1..8],'Gene_pheno_score','Mix_score',@head[9..$#head]);
#@out_head=("#$head[0]",@head[1..6],@head[7..$#head]);

##main shell

#set path to outdir if exists outdir
if($outdir){
    mkpath $outdir;
    $path=abs_path($outdir);
}

##get each variant info,get all samples and each sample genotype info;
readAnnoVcfList(\%hash,$vcf) if $vcf;
readAllInfo($vcflist) if $vcflist;

##Add other database;
readGeneScore($genescore);
readNSFP_Gene(\%hNSFP_Gene,$dbNSFP_gene);

#interp file
my($interpFILE,$revealDB,$mcapDB,$dbNSFPDB,$dbSNPDB);
if($interp){
    die("$interp suffix should be gz or list") if $interp!~/gz|list$/i;
    if($interp=~/list/i){
        open IN,$interp or die "not open $interp file\n";
        open OUT,">$path/interp.tsv" or die $!;
        my %hash_interp;
        my @interp_list=<IN>;
        for(my $i=0;$i<=$#interp_list;$i++){
            if($interp_list[$i]=~/.gz$/i){
                open FILE,"gzip -dc $interp_list[$i]|" or die $!;
            }else{
                open FILE,"$interp_list[$i]" or die $!;
            }

            while(<FILE>){
                chomp;
                my @info=split;
                if(not exists $hash_interp{$info[0]}{$info[1]}{$info[3].$info[4]}){
                    $hash_interp{$info[0]}{$info[1]}{$info[3].$info[4]}[0]=1;
                    $hash_interp{$info[0]}{$info[1]}{$info[3].$info[4]}[1]=join "\t",@info;
                }
            }   
        }
        
        
         for my $chr(1..22,'X','Y'){
           foreach my $pos (sort{$a<=>$b} keys $hash_interp{$chr}){
                foreach my $gt(keys $hash_interp{$chr}{$pos}){
                   print OUT "$hash_interp{$chr}{$pos}{$gt}[1]\n";          
                }
            }
        }
        
        system "$bgzip -f $path/interp.tsv && $tabix -s 1 -b 2 -e 2 $path/interp.tsv.gz";
        $interpFILE  = Tabix->new(-data=>"$path/interp.tsv.gz")
    }else{
        $interpFILE  = Tabix->new(-data=>$interp);
    }
}
$revealDB    = Tabix->new(-data=>$reavel);
$mcapDB      = Tabix->new(-data=>$mcap);
$dbSNPDB    = Tabix->new(-data=>$dbSNP) if ($dbsnp);

AddDatabase(\%hash);

#group and print
if($pedgree){
    readPedInfo($pedgree);
    
    getAllSamplesInfo(\%hash,\@AllSample_vcf);

    foreach my $eachfam (keys %fam_n){
        printData(\%hash,\%hNSFP_Gene,$eachfam);
    }
}else{
   #if heven't pedgree,set all sample to case;
   for((@AllSample_vcf)){
       $fam_n{'all'}{$_}=2;
   }
   #@{$sample_fam{'all'}}=@AllSample_vcf; 
   getAllSamplesInfo(\%hash,\@AllSample_vcf);
   printData(\%hash,\%hNSFP_Gene,'all');
}

#system "iconv -c -f utf-8 -t gbk $path/$name.result_explain_filtered.tsv >$path/$name.result_explain_filtered.win.xls";
getStat(\%hash,'SNP','1');#outfile,varant type,pop freq
printStat("$path/$sample.Anno_stat_snp_all.xls",'SNP');
getStat(\%hash,'SNP','0.05');#outfile,varant type,pop freq
printStat("$path/$sample.Anno_stat_snp_frqLt0.05.xls",'SNP');
getStat(\%hash,'InDel','1');#outfile,varant type,pop freq
printStat("$path/$sample.Anno_stat_indel_all.xls",'InDel');
getStat(\%hash,'InDel','0.05');#outfile,varant type,pop freq
printStat("$path/$sample.Anno_stat_indel_frqLt0.05.xls",'InDel');

##readAllInfo
sub readAllInfo{
    my ($vcflist)=shift;
    open IN,"$vcflist"or die $!;
    my @info=<IN>;
    for(my $i=0;$i<=$#info;$i++){
        my $vcf=(split(/\t/,$info[$i]))[1];    
        readAnnoVcfList(\%hash,$vcf);
    }
}


##readInfo
sub readPedInfo{
    my ($pedigree)=shift;
    
    open IN,"$pedigree" or die $!;
    while(<IN>){
        chomp;
        next if /^#/;
        next if /^$/;
        my @info_line=split;
        my($fam,$indv,$phe)=@info_line[0,1,5];
        
        $fam_n{$fam}{$indv}=$phe;
     }
}


sub readDieseaInfo{
    my($dieaseinfo)=@_;
    open IN,$dieaseinfo or die $!;
    while(<IN>){
        chomp;
        my @line=split(/\t/,$_);
        @{$dieaseinfo{$line[0]}}=@line[1..$#line];
    }
    
    print &get_time(10)."\tDiease information file $dieaseinfo read finished\n";
}

sub readCardInfo{
    my($Cardinfo)=@_;
    if($Cardinfo=~/gz$/){
        open IN,"gunzip -dc $Cardinfo|" or die $!;
    }else{
        open IN,"$Cardinfo" or die $!;
    }
    while(<IN>){
        chomp;
        my @line=split(/\t/,$_);
        @{$Cardinfo{$line[0].$line[1].$line[3].$line[4]}}=@line[5..$#line];
    }
    print &get_time(10)."\tCard information $Cardinfo read finished\n";
}


sub readGeneScore{
     my($genescore)=@_;
     if($genescore=~/gz$/){
        open IN,"gunzip -dc $genescore|" or die $!;
     }else{
        open IN,"$genescore" or die $!;
     }
     while(<IN>){
        chomp;
        my @line=split(/\t/,$_);
        $genescore{$line[1]}=$line[3]
    }
    print &get_time(10)."\tDisease Gene Score $genescore read finished\n";

}

sub readNSFP_Gene{
    my($hNSFP_Gene,$f_nfsp)=@_;
    
    open NFSP,"gzip -dc $f_nfsp|" or die $!;
    <NFSP>;
    while(<NFSP>){
        chomp;
        next if /^$/;
        my @info=split(/\t/,$_);
        @{$$hNSFP_Gene{$info[0]}} = (@info[17..20,24..26]); #gene=>(Pathway(KEGG)_id,Pathway(KEGG)_full,Function_description,Disease_description,GO_biological_process,GO_cellular_component,GO_molecular_function)
    }
    print &get_time(10)."\tNSFP_Gene $f_nfsp read finished\n";
}


sub readAnnoVcfList{
    my($hash,$anno_vcf)=@_;

    #open file
    if($anno_vcf=~/gz$/){
        open IN,"gzip -dc $anno_vcf|" or die "not open anno vcf file $anno_vcf\n";
    }else{
        open IN,$anno_vcf or die "not open anno vcf file $anno_vcf\n";
    }

    #get prefix and path 
    my @suffixlist=qw/.vcf .vcf.gz/;
    ($name,$path)=fileparse($anno_vcf,@suffixlist);

    #get sample name;

   if(!$sample){
       $sample=(split(/\./,$name))[0];
   }
   
   my @sample_eachVCF;
   while(<IN>){
       chomp;
       if($_=~/^#CH/){
           my @head_info=split(/\t/,$_);
           push @AllSample_vcf,@head_info[9..$#head_info];
           push @sample_eachVCF,@head_info[9..$#head_info];
       }
       
       next if /^#/;
       my $info=$_;
       my @line=split(/\t/,$info);
       my ($chr, $start, $ID, $ref_allele, $mut_allele, $quality_score, $filter, $info_detail, $format, @sample_gt) = @line;
       my @format=split(/:/,$line[8]);
       my @format_index;
       
       my ($LOF,$NMD,$zygosity,$alt_depth,$exac_eas_af,$kg_eas_af,$genotype,$variant_type);

#snp and indel type
       next if ($mut_allele eq '.');
   

        my $mut_len;
        if($mut_allele=~/,/){
           $mut_len = length ((split(/,/,$mut_allele))[0]);
        }else{
           $mut_len = length($mut_allele);
        }
 
        if (length($ref_allele)==1 && $mut_len ==1 && $mut_allele =~m/[ATCGatcg]/){
            $variant_type="SNP";
        } else {
            $variant_type="InDel";
        }


#handle ANN info        
#Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
##C|missense_variant|MODERATE|F5|ENSG00000198734|transcript|ENST00000367797|protein_coding|13/25|c.3980A>G|p.His1327Arg|4182/7024|3980/6675|1327/2224||;
       
           my ($gene,$high_effect,$high_effect_ID,$impact,$multi_gene,$annotation,$Featrre_ID,$rank_exon,$HGVS_c,$HGVS_p)=('.','.','.','.','.','.','.','.','.','.');
           if ($info=~/ANN=([^\;]+)/){
               my $anno=$1;
               ($gene,$high_effect,$high_effect_ID,$impact,$multi_gene,$annotation,$Featrre_ID,$rank_exon,$HGVS_c,$HGVS_p)=FetchAnnoInfo($anno);
           }
           else{
               next
           }
              
#handle LOF and NMD      
           ($info=~/LOF=[^\|]+\|[^\|]+\|([^\)]+)/) ? ($LOF=$1) : ($LOF=$na_mark);#LOF=(OBSCN|ENSG00000154358|3|1.00) 
           ($info=~/NMD=[^\|]+\|[^\|]+\|([^\)]+)/) ? ($NMD=$1) : ($NMD=$na_mark);#NMD=(OBSCN|ENSG00000154358|3|1.00)

#handle freq
           ($info=~/dbNSFP_ExAC_EAS_AF=([^;]+)/) ?($exac_eas_af=$1):($exac_eas_af=$na_mark);
           ($info=~/dbNSFP_1000Gp3_EAS_AF=([^;]+)/) ?($kg_eas_af=$1):($kg_eas_af=$na_mark);
           
           $chr=~s/chr//;
   
           my $titv=titv($ref_allele, $mut_allele);
     

      
         @{$$hash{$chr}{$start}{$ref_allele}{$mut_allele}{'titv'}} = $titv;
         @{$$hash{$chr}{$start}{$ref_allele}{$mut_allele}{'anno'}} = ($multi_gene,$annotation,$Featrre_ID,$gene,$high_effect,$high_effect_ID,$HGVS_p,$rank_exon,$HGVS_c,$impact,$LOF,$NMD);#anno info; chr_pos_ref_mut =>(gene,all anno function,all anno transcript ID ,high effect function,high effect transcript ID, HGVS_p,HGVS_c,Imapct mark,lof,nmd )  
         @{$$hash{$chr}{$start}{$ref_allele}{$mut_allele}{'id'}}   = ($ID,$variant_type);
         

         for (my $i=0;$i<=$#sample_eachVCF;$i++){
                 push @{$$hash{$chr}{$start}{$ref_allele}{$mut_allele}{$sample_eachVCF[$i]}},$sample_gt[$i];
         }
     }
         print &get_time(10)."\tanno vcf inforamtion file $anno_vcf read finished\n";
}


sub getAllSamplesInfo{
    my ($hash,$AllSample_vcf)=@_;
    my @AllSample_vcf=@$AllSample_vcf;
    my $flag=1;
    #for my $chr(1..22, 'X', 'Y', 'MT'){ 
     for my $chr(sort{$a cmp $b} keys %$hash){ 
       for my $pos ( sort{$a <=> $b} keys %{$$hash{$chr}}){
            for my $ref_allele (keys %{$$hash{$chr}{$pos}}){
                for my $mut_allele (keys %{$$hash{$chr}{$pos}{$ref_allele}}){
                    #trans pedgree info to array position  
                    foreach my $fam( sort{$a cmp $b} keys %fam_n){
                        #my @sample_fam; 
                        my @sample_eachfam_gt;
                         
                        for(my$i=0;$i<=$#AllSample_vcf;$i++){
                            if(exists $fam_n{$fam}{$AllSample_vcf[$i]}){
                                push @{$sample_fam{$fam}},$AllSample_vcf[$i] if $flag==1;
                                if(exists $$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$AllSample_vcf[$i]}){
                                    push @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'gt'}},
                                         @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$AllSample_vcf[$i]}};
                                    push @sample_eachfam_gt,
                                         @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$AllSample_vcf[$i]}};
                                }else{
                                    push @sample_eachfam_gt,'0/0';
                                }
                            }
                        }
                       
                        #get family info
                                                   

                            my ($zygosity,$all_alt,$alles_info,$alt_freq,$twotailed_value)
                                =SamZygJuge(\@sample_eachfam_gt,$fam,\@{$sample_fam{$fam}},$ref_allele,$mut_allele);

                            @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'zyg'}} 
                                = ($zygosity,$all_alt,$alles_info,$alt_freq,$twotailed_value);
                            
                            my ($mut_sample,$mut_info)
                                =getmMutSample(\@sample_eachfam_gt,$fam,\@{$sample_fam{$fam}},$ref_allele,$mut_allele);

                            @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'mutsample'}}
                                =($mut_sample,$mut_info);  
                            @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'gt'}}=@sample_eachfam_gt;
                            push @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'allgt'}},@sample_eachfam_gt;
                    }
                        $flag=0;
                 
                }
            }
        }
    }
}

sub FetchAnnoInfo{
    my ($anno) =@_;
    my ($multi_gene,$high_effect,$impact,$gene,$annotation,$Featrre_ID,$rank_exon,$HGVS_c,$HGVS_p,$high_effect_ID);
    
    $anno=(split(/\t/,$anno))[0];
    my @anno=split(/,/,$anno);
    my $i=0;
    for((@anno)){
        my @anno_info=split(/\|/,$_);
        $i++;
        
        
        my @array_effect=split(/&/,$anno_info[1]);
      
        my $effect = $array_effect[0];
            
         
        for(my $j=0;$j<=$#array_effect;$j++){

            if(not exists $RankScore{$array_effect[$j]}){
                if($anno_info[2]=~/High/i){
                    $RankScore{$array_effect[$j]}=8;
                }elsif($anno_info[2]=~/MODERATE/i){
                    $RankScore{$array_effect[$j]}=5;
                }elsif($anno_info[2]=~/MODIFIER/i){
                    $RankScore{$array_effect[$j]}=4;
                }elsif($anno_info[2]=~/LOW/i){
                    $RankScore{$array_effect[$j]}=2;
                }
            }
           
           $effect=$array_effect[$j] if $RankScore{$array_effect[$j]}>$RankScore{$effect};
        }

 ##get high effect transcrpit           

 
            if($i==1){
                $annotation     = $anno_info[1];
                $Featrre_ID     = $anno_info[6];
                $multi_gene     = $anno_info[3];
                $high_effect    = $effect;
                $high_effect_ID = $anno_info[6];
                $rank_exon      = judgeExist($anno_info[8]);
                $HGVS_c         = judgeExist($anno_info[9]);
                $HGVS_p         = judgeExist($anno_info[10]);
                $impact         = $anno_info[2];
                $gene           = judgeExist($anno_info[3]);
 
            }else{
                if($RankScore{$effect} gt $RankScore{$high_effect}){
                    $high_effect    = $effect;
                    $high_effect_ID = $anno_info[6];
                    $impact         = $anno_info[2];
                    $rank_exon      = judgeExist($anno_info[8]);
                    $HGVS_c         = judgeExist($anno_info[9]);
                    $HGVS_p         = judgeExist($anno_info[10]);
                    $gene           = judgeExist($anno_info[3]);
                }
               $annotation.=",$anno_info[1]";
               $Featrre_ID.=",$anno_info[6]";
               $multi_gene .=",$anno_info[3]"; 
           }

        }
      return ($gene,$high_effect,$high_effect_ID,$impact,$multi_gene,$annotation,$Featrre_ID,$rank_exon,$HGVS_c,$HGVS_p);
}



sub AddDatabase{
    my ($hash,$type,$pop_freq)=@_;
    for my $chr(1..22, 'X', 'Y', 'MT'){ 
       for my $pos ( sort{$a <=> $b} keys %{$$hash{$chr}}){
            for my $ref_allele (keys %{$$hash{$chr}{$pos}}){
                for my $mut_allele (keys %{$$hash{$chr}{$pos}{$ref_allele}}){
                    my ($gene,$annotation,$Featrre_ID,$high_effect,$high_effect_ID,$HGVS_p,$HGVS_c,$impact,$LOF,$NMD)
                        =@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'anno'}}; 
                        #anno info; chr_pos_ref_mut =>(gene,all anno function,all anno transcript ID ,
                        #high effect function,high effect transcript ID, HGVS_p,HGVS_c,Imapct mark,lof,nmd )  

#MCAP info
                    my ($mcap_score,$mcap_pred)=(".",".");
                        $mcap_score=getMCAP($chr,$pos-1,$pos,$ref_allele, $mut_allele);
             
                    if($mcap_score=~/^\./){
                        $mcap_pred=".";
                    }elsif($mcap_score>=0.025){
                        $mcap_pred="Possible pathogenic";
                    }elsif($mcap_score<0.025){
                        $mcap_pred="Possible benign";
                    }else{
                        $mcap_pred="uncertain";
                    }
        
#REOHIT result
                   my ($interp_content);
                   ($interp_content)=getREOHIT($chr,$pos-1,$pos,$ref_allele, $mut_allele) if ($interp);

##Add REAVEL

                   my ($reaveal_score)=(".");   
                   my $reaveal_pred ;        
          
                   $reaveal_score=getREVEL($chr,$pos-1,$pos,$ref_allele, $mut_allele);

                   if($reaveal_score=~/^\./){
                       $reaveal_pred=".";
                   }elsif($reaveal_score>=0.95){
                       $reaveal_pred="Possible pathogenic";
                   }elsif($reaveal_score<=0.05){
                       $reaveal_pred="Possible benign";
                   }else{
                       $reaveal_pred="uncertain";
                   }

##Add dbsnp
                 my $rs=(".");
                 if($dbsnp){
                     $rs = getDBsnp($chr,$pos-1,$pos,$ref_allele, $mut_allele);
                 }
                 

#count genescore and impact score

#genesore
                  my $genescore=".";
                  if(exists $genescore{$gene}){
                      $genescore=$genescore{$gene};
                  }else{
                      $genescore=0;
                  }

##impact sore

                    my $impactscore;
                    my($freqscore,$funcscore,$predscore,$thirdPartScore,$totalscore);
                    
                    #1000g Allele Frequency(EAS),ExAC Allele Frequency(EAS),In-house Allele Frequency(EAS);
                    $freqscore = popFreqScore($$interp_content[1],$$interp_content[4],$$interp_content[5]);
                    
                    $funcscore = funcScore($RankScore{$high_effect},$LOF,$NMD);
                   
                    #Missense prediction result,Affect splicing,Highly conserved,
                    $predscore = predScore($$interp_content[23],$$interp_content[24],$$interp_content[25]);
                    
                    $thirdPartScore =thirdPartScore($mcap_score,$reaveal_score);
                    
                    $totalscore=$freqscore+$funcscore+$predscore+$thirdPartScore;
                    $max_impact_score=$totalscore if $totalscore>$max_impact_score;
                  


        #push Data and score;
        @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'db'}}=(@$interp_content[0..5],$mcap_score,$mcap_pred,$reaveal_score,$reaveal_pred,@$interp_content[23..25],@$interp_content[6..22],@$interp_content[26..$#$interp_content]);
        @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'score'}}=($totalscore,$genescore);
        #$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'rs'}=$rs;
        @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'id'}}[0] = $rs if $dbsnp;

        #@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}[1]}=(@$interp_content[0..5],$mcap_score,$mcap_pred,$reaveal_score,$reaveal_pred,@$interp_content[22..24],@$interp_content[6..21],@$interp_content[25..$#$interp_content],$freqscore,$funcscore,$predscore,$thirdPartScore,$totalscore);
        # print join ("\t",$chr,$pos,$ref_allele, $mut_allele,$freqscore,$$interp_content[1],$$interp_content[4],$$interp_content[5]),"\n";
        # print join ("\t",$high_effect,$RankScore{$high_effect},$LOF,$NMD,$funcscore),"\n";
        # print join ("\t",$predscore,$$interp_content[22],$$interp_content[23],$$interp_content[24]),"\n";
        # print join ("\t",$thirdPartScore,$mcap_score,$mcap_pred,$reaveal_score,$reaveal_pred),"\n";
         }
       }
     }
  }
        print &get_time(10)."\tAll data added and score evaluate finished\n";
}


sub popFreqScore{
    my($kg_eas,$exac_eas,$inhouse)=@_;
    
     
    my $score=0;
    my $record=0;
  
   #0.003	0.00323624595469	. 
    
   for(($kg_eas,$exac_eas)){
        $_=~s/^\./0/;
        if($_ gt 0.05){#commom variant;
            $record-=1;
        }elsif($_ lt 0.0001){#less than 0.01%,vary rare
            $record+=1.25; 
        }elsif($_ lt 0.005){#less than 0.5%(1/500)
            $record+=1;
        }elsif($_ lt 0.01){#less than 1%,rare variant
            $record+=0.5;
        }elsif($_ lt 0.05){#less than 5%,rare variant
            $record+=0.25;
        }
     }
   
 
    if($inhouse){
        for(($inhouse)){
            $_=~s/^\./0/;
            if($_ gt 0.15){#very common in house database,15%;minus score
                $record-=0.5; 
            }elsif($_ gt 0.01 && $_ lt 0.15){#common in house database,1%-15%,minus score,eg,0.15*-1
                $record = $record-$_; 
            }
         }
     }
 
   #turn record to score;
    if($record>1){
        $score=0.75+($record-1)*0.25;
    }elsif($record<-1){
        $score=-0.75+($record+1)*0.25;
    }else{
        $score=$record; 
    }
    #turn to standar_score
    $score=$standard_score{'PS'}*$score;
    return ($score); 
}

sub funcScore{
    my ($func_raw_score,$LOF,$NMD)=@_;#func raw score,$lof info,$nmd info
    my $score=0;
    if ($NMD!~/^\./ && $LOF!~/^\./){
        $func_raw_score+=1.5;
    }elsif ($NMD!~/^\./ || $LOF!~/^\./){
        $func_raw_score+=1;
    }

    $score = $func_raw_score/12*$standard_score{'PVS'};
    return($score);
} 

sub predScore{
    my($missense_pred,$splice_pred,$conserv_pred)=@_;
    my ($score,$impactscore)=(0,0);
  
#    13.05   Damaging        Yes     Yes
 
    $impactscore+=1 if $missense_pred =~ /Damaging/;
    $impactscore-=0.5 if $missense_pred =~ /Tolerated/;
     
    
    $impactscore+=1 if $splice_pred =~ /Yes/;

    $impactscore+=1 if $conserv_pred =~ /Yes/;
    $impactscore-=0.5 if $conserv_pred =~ /NO/;
   
    
    if($impactscore>=2){
        $score=1+($impactscore-1)*0.25;
    }elsif($impactscore<=-2){
        $score=-1+($impactscore+1)*0.25;
    }else{
        $score=$impactscore;
    }
 
    $score=$standard_score{'PP'}*$score;
    return ($score);    
}

sub thirdPartScore{
    my($mcap_score,$reaveal_score)=@_;
    
    my ($score,$impactscore)=(0,0);
    

    #high reaveal_sore,high imapct.4levels
    if($reaveal_score=~/^\./){
        $impactscore+=0;
    }elsif($reaveal_score gt 0.95){
        $impactscore+=$reaveal_score;
    }elsif($reaveal_score lt 0.05){
        $impactscore=$impactscore-$reaveal_score-1;
    }else{
        $impactscore+=$reaveal_score/2;
    }   
   
   #exists mcap_score and mcap_score less than 0.0125 
    if($mcap_score =~/^\./){
        $impactscore+=0;
    }elsif($mcap_score le 0.0125){
        $impactscore=$impactscore-$mcap_score-1;
    }else{
        $impactscore=$impactscore+$mcap_score+0.5;
    }

    if($impactscore>1){
        $score=1+($impactscore-1)*0.5;
    }elsif($impactscore<-1){
        $score=-1+($impactscore+1)*0.5;
    }else{
        $score=$impactscore;    
    }
    
    $score=$standard_score{'PP'}*$score;
    return ($score);
}


sub printData{
    my ($hash,$hNSFP_Gene,$fam)=@_;
 
    my $out_fam=""; 
       $out_fam="_$fam" if $fam!~/all/;

    if($vcflist){ 
       $name="combine";
    }
    open OUT,">$path/$name.result_explain$out_fam.tsv" or die $!;
    open OUT1,">$path/$name.result_explain_filtered$out_fam.tsv" or die $!;
        
    print OUT join ("\t",@out_head,@{$sample_fam{$fam}}),"\n";
    print OUT1 join ("\t",@out_head,@{$sample_fam{$fam}}),"\n";
     
    for my $chr(1..22, 'X', 'Y', 'MT'){ 
       for my $pos ( sort{$a <=> $b} keys %{$$hash{$chr}}){
            for my $ref_allele (keys %{$$hash{$chr}{$pos}}){
                for my $mut_allele (keys %{$$hash{$chr}{$pos}{$ref_allele}}){
                   my ($mut_sample,$mut_info) = @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'mutsample'}};

                   my ($impactscore,$genescore)=@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'score'}};
                   my $normalized_impac_score  = $impactscore/$max_impact_score;
                   my $mixscore = $normalized_impac_score*0.65+$genescore*0.35;
                   
                   my ($multi_gene,$annotation,$Featrre_ID,$gene,$high_effect,$high_effect_ID,$HGVS_p,$rank_exon,$HGVS_c,$impact,$LOF,$NMD)
                       = @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'anno'}};
                   #anno info; chr_pos_ref_mut =>(gene,all anno function,all anno transcript ID ,high effect function,
                   #high effect transcript ID, HGVS_p,HGVS_c,Imapct mark,lof,nmd )  
                   my ($ID,$variantType)=@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'id'}};
                   #$ID=$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'rs'} if $dbsnp;

                   my ($kg_eas_af,$exac_eas_af,$in_house_af) 
                      = @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'db'}}[1,4,5];
                   #1000g Allele Frequency(EAS),ExAC Allele Frequency(EAS),In-house Allele Frequency(EAS);

                   my ($zyg) = $$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'zyg'};
                   #add gene GO and KO
                   my @gene_function;
                   if(exists $$hNSFP_Gene{$gene}){
                       @gene_function=@{$$hNSFP_Gene{$gene}};
                   }else{
                       @gene_function=qw/. . . . . . ./;
                   } 
      
                   next if $zyg=~/NONEMUT/;             
                   print  OUT join ("\t",
                    $mut_sample,$mut_info,
                    $chr,$pos,$ID,$ref_allele,$mut_allele,$variantType,$normalized_impac_score,$genescore,$mixscore,
                    @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'zyg'}},
                    @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'anno'}},
                    @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'db'}},
                    @gene_function,@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'gt'}}
                    ),"\n";        
                  
                  #repalce allele frequency '.' to '0'
                  for($kg_eas_af,$exac_eas_af,$in_house_af){
                     if($_=~/^.$/){
                        $_=~s/./0/;
                     }
                  }                  

                  ## write filter file
                  if($RankScore{$high_effect} > 5 && 
                     $impactscore > 55 &&
                     $kg_eas_af < $kg_eas_af_min &&
                     $exac_eas_af < $exac_eas_af_min &&
                     $in_house_af < $in_house_af_min
                     )
                                      
                  {
                      print  OUT1 join ("\t",
                        $mut_sample,$mut_info,                   
                        $chr,$pos,$ID,$ref_allele,$mut_allele,$variantType,$normalized_impac_score,$genescore,$mixscore,
                        @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'zyg'}},
                        @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'anno'}},
                        @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'db'}},
                        @gene_function,@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{$fam}{'gt'}}
                        ),"\n";        
                  }
              }
          }
      }
  } 
    close OUT;
}
 

sub getStat{
    my ($hash,$type,$pop_freq)=@_;
    #my(@{$$hash{$chr.$start.$ref_allele.$mut_allele}[0]})=@_;
    ($nstopg,$nstopl,$nstartl,$nsplicea,$nspliced,$nframe,$ndbsnp,$nnovel,$inframe_ins,$inframe_del,$inframe_dirup_ins,
      $inframe_dirup_del,$nintergenic,$nexonice,$nsplicing,$nUTR,$nsplice_region,$nintronic,$nstream,
      $nnoncoding,$nTFbind,$nintragenic,$ti,$tv,$nvariants,$nmissense,$nsynonymous,$nlof,$nnmd)=
    (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

    for my $chr(1..22, 'X', 'Y', 'MT'){ 
       for my $pos ( sort{$a <=> $b} keys %{$$hash{$chr}}){
            for my $ref_allele (keys %{$$hash{$chr}{$pos}}){
                for my $mut_allele (keys %{$$hash{$chr}{$pos}{$ref_allele}}){
                    my ($ID,$variant_type)=@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'id'}};
                    my ($gene,$annotation,$Featrre_ID,$high_effect,$high_effect_ID,$HGVS_p,$rank_exon,$HGVS_c,$impact,$LOF,$NMD)=@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'anno'}};
                    my($titv) = @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'titv'}};
                    
                    my($kg_eas_af) = "."; 
                    $kg_eas_af = @{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'db'}}[1] if(@{$$hash{$chr}{$pos}{$ref_allele}{$mut_allele}{'db'}}[1]);

                    $kg_eas_af=(split(/,/,$kg_eas_af))[0] if $kg_eas_af=~/\,/;
                    
                    if($type=~/$variant_type/i){
                        if($kg_eas_af=~/^\./ || $kg_eas_af lt $pop_freq){
                            funcStat($annotation,$HGVS_p);
                            $nlof++ if $LOF!~/^\./;
                            $nnmd++ if $NMD!~/^\./;
                            if($ID!~/^\./){
                                $ndbsnp++;
                            }else{
                                $nnovel++;
                            }

                            if($type=~/snp/i){
                                $ti++ if $titv=~/ti/i;
                                $tv++ if $titv=~/tv/i;
                           }                  
                       }
                   }
                }
            }   
       }
    }
}

sub getREVEL{
    my($chr,$start,$end,$ref_allele,$mut_allele)=@_;

    my ($reaveal_score)=(".");   
        
    if($chr=~/^[0-9]+|^[XY]/){
    my $revealIDX = $revealDB->query($chr,$start-1,$start);

    while(my $record=$revealDB->read($revealIDX)){
         
        my @reveal_line=split(/\t/,$record);

        if ($chr eq "$reveal_line[0]" && $start eq "$reveal_line[1]" && $ref_allele eq $reveal_line[2] && $mut_allele eq $reveal_line[3] ){
             $reaveal_score=$reveal_line[6];
         }
    }
    }
    return $reaveal_score;
}

sub getMCAP{
    my($chr,$start,$end,$ref_allele,$mut_allele)=@_;

    my ($mcap_score)=(".");   
      
    if($chr=~/^[0-9]+|^[XY]/){
        
        my $mcapIDX = $mcapDB->query($chr,$start,$end);
        
        while(my $record=$mcapDB->read($mcapIDX)){

            my @mcap_line=split(/\t/,$record);

            if ($chr eq "$mcap_line[0]" && $end eq "$mcap_line[1]" && $ref_allele eq $mcap_line[2] && $mut_allele eq $mcap_line[3] ){
                $mcap_score=$mcap_line[4];
            }
        }
    }
    return $mcap_score;
    
}


sub getREOHIT{
    my($chr,$start,$end,$ref_allele,$mut_allele)=@_;

    my @interp_result;

    if($chr=~/[0-9]+|[XYM]/){
        my $interpIDX = $interpFILE->query($chr,$start,$end);

        while(my $record=$interpFILE->read($interpIDX)){
 
            my @line=split(/\t/,$record);
          
            if ($chr eq "$line[0]" && $end eq "$line[1]" && $ref_allele eq $line[3] && $mut_allele eq $line[4]){
                push @interp_result,@line[19,20,22,25,26,28..45,58,61,68..73];
                #1000g Allele Frequency,1000g Allele Frequency(EAS),EVS Allele Frequency,ExAC Allele Frequency,ExAC Allele Frequency(EAS),In-house Allele Frequency(EAS);1-6->19,20,22,25,26,28.
                #CGD Condition(gene),CGD Inheritance,CGD PumMed,OMIM Disease(gene),OMIM Gene,ClinVar_ID,ClinVar Clinical Significance,ClinVar Report,ClinVar Star,ClinVar Disease,HGMD ID,HGMD Categorization,HGMD Disease,HGMD PubMed,GWAS OR/Beta,GWAS trait,GWAS PubMed;7-22 -> 29-45
                #Missense prediction result,Affect splicing,Highly conserved,Repeat name,Functional domain name,Category,evidences,evidences information;24-31 ->58,61,68-73
            }
        }
    }
    return (\@interp_result);
}


sub getDBsnp{
    my($chr,$start,$end,$ref_allele,$mut_allele)=@_;

    my $rs=".";
    if($chr=~/[0-9]+|[XYM]/){
       $chr="chr$chr";#for vcf4.0
       $chr="chrM" if $chr=~/MT/; 
       my $dbsnpIDX = $dbSNPDB->query($chr,$start,$end);

        while(my $record=$dbSNPDB->read($dbsnpIDX)){

            my @line=split(/\t/,$record);
            if ($chr eq "$line[0]" && $end eq "$line[1]" && $ref_allele eq $line[3] && $line[4]=~/$mut_allele/){
                $rs=$line[2];
            }
        }
    }
    return ($rs);
}


sub judgeExist{
    my($input)=@_;
    if($input){
        if($input=~/\S+/){
            return $input;
        }else{
            return $na_mark;
        }
    }else{
        return $na_mark;
    }
}

sub funcStat{
    my ($func,$HGVS_p)=@_;
    $nvariants++;
#stop gain and stop loss
    if ($func =~ /\bstop_gained\b/){$nstopg++;}
    if ($func =~ /\bstop_lost\b/){$nstopl++;}
    if ($func =~ /\bstart_lost\b/){$nstartl++;}
    
    if ($func =~ /\bmissense/){$nmissense++;}
    elsif ($func =~ /\bsynonymous/){$nsynonymous++;}


##splice
       if($func=~/\bsplice_acceptor_variant\b/){$nsplicea++};
       if($func=~/\bsplice_donor_variant\b/){$nspliced++};
       
#frame
     if ($func =~ /\bframeshift_variant\b/){$nframe++;}
     if ($func =~ /\binframe_insertion\b/){$inframe_ins++;}
     if ($func =~ /\binframe_deletion\b/){$inframe_del++;}
     if ($func =~ /\bdisruptive_inframe_deletion\b/){$inframe_dirup_del++;}
     if ($func =~ /\bdisruptive_inframe_insertion\b/){$inframe_dirup_ins++;}

#region func
     if($HGVS_p=~/^p/){$nexonice++;}
     elsif($func=~/\bprotein_protein_contact\b/){$nexonice++;}
     elsif($func=~/\bintergenic/){$nintergenic++;}
     elsif($func=~/\b(splice_donor_variant|splice_acceptor_variant)\b/){$nsplicing++;}
     elsif($func=~/\bintron/){$nintronic++;}
     elsif($func=~/prime_UTR|initiator_codon_variant/){$nUTR++;}
     elsif($func=~/\b(upstream|downstream)/){$nstream++;}
     elsif($func=~/\b(splice_region_variant)\b/){$nsplice_region++;}
     elsif($func=~/\b(non_coding_exon_variant)\b/){$nnoncoding++;}
     elsif($func=~/\b(intragenic_variant)\b/){$nintragenic++;}
     elsif($func=~/\b(TF_binding_site_variant)\b/){$nTFbind++;}
}

sub titv{
    my($base1,$base2)=@_;
    if((($base1 eq 'A') && ($base2 eq 'G')) || (($base1 eq 'G') && ($base2 eq 'A')) || (($base1 eq 'C') && ($base2 eq 'T')) || (($base1 eq 'T') && ($base2 eq 'C'))){
      return "ti";
    }else{
      return "tv";
    }
}

sub printStat{
    my($outfile,$vtype)=@_;
    my($ti_tv);

    open OUT,">$outfile" or die $!;
    print OUT "Sample\t$sample\n";
    print OUT "Variant_type\t$vtype\n";
    print OUT "total_Variants\t$nvariants\n";
    printf OUT ("Dbsnp\t%d(%0.2f%%)\n",$ndbsnp,$ndbsnp/$nvariants*100);
    printf OUT ("Noval_snps\t%d(%0.2f%%)\n",$nnovel,$nnovel/$nvariants*100);
    print OUT "Stop_gain\t$nstopg\n";
    print OUT "Stop_loss\t$nstopl\n";
    print OUT "Start_lost\t$nstartl\n";
    print OUT "LOF(loss_of_function)\t$nlof\n";
    print OUT "NMD(nonsense_mediated_decay)\t$nnmd\n";
    
    if($vtype=~/snp/i){
        $ti_tv=sprintf("%.2f",$ti/$tv);
        print OUT "Missense\t$nmissense\n";
        print OUT "Synonymous\t$nsynonymous\n";
    }
    
    if($vtype=~/indel/i){
        print OUT "frameshift\t$nframe\n";
        print OUT "inframe_deletion\t$inframe_del\n";
        print OUT "inframe_insertion\t$inframe_ins\n";
        print OUT "disruptive_inframe_deletion\t$inframe_dirup_del\n";
        print OUT "disruptive_inframe_insertion\t$inframe_dirup_ins\n";
    }
    
    print OUT "Exonic\t$nexonice\n";
    print OUT "Splicing_acceptor\t$nsplicea\n";
    print OUT "Splicing_donor\t$nspliced\n";
    print OUT "Splicing_region\t$nsplice_region\n";
    print OUT "Non_coding_exon\t$nnoncoding\n";
    print OUT "UTR\t$nUTR\n";
    print OUT "Intronic\t$nintronic\n";
    print OUT "Upstream&downstreams\t$nstream\n";
    print OUT "Intergenic\t$nintergenic\n";
    print OUT "Intragenic\t$nintragenic\n";
    print OUT "TF_binding_site\t$nTFbind\n" if $vtype=~/snp/i;
    print OUT "Ratio_Ti/Tv\t$ti_tv\n" if $vtype=~/snp/i;
}

sub get_time {
    my $interval = $_[0]*60;

    my ($sec,$min,$hour,$day,$mon,$year,$weekday,$yeardate,$savinglightday) = (localtime(time + $interval));

    $sec  = ($sec < 10)? "0$sec":$sec;
    $min  = ($min < 10)? "0$min":$min;
    $hour = ($hour < 10)? "0$hour":$hour;
    $day  = ($day < 10)? "0$day":$day;
    $mon  = ($mon < 9)? "0".($mon+1):($mon+1);
    $year += 1900;
  
    return "$year-$mon-$day $hour:$min:$sec.00";
}


sub Zygosity{
    my ($info,$ref,$alt,$format)=@_;
     
    my (@format_index,$zygosity,$genotype,$alt_depth);
    if($format){
        my @format=split(/:/,$format);
        for(my $i=0;$i<=$#format;$i++){
            push @format_index,$i if $format[$i]=~/\b(GT|AD|DV|GQ|PL)\b/;
        }
    }else{
        @format_index=qw/0 1 2/;
    }

    #my ($GT,$AD,$GQ)=(split(/:/,$info))[@format_index];
    my ($GT,$AD)=(split(/:/,$info))[0,1];

    my @alt;
    if($alt=~/,/){
        @alt=split(/,/,$alt);
    }else{
        push @alt,$alt;
    }

    if ($GT =~ m#^0/1# || $GT =~ m#^1/0#){
        $zygosity = "het";
        if($AD=~/,/){
            $alt_depth=(split(/,/,$AD))[1];
        }else{
            $alt_depth=$AD;
        }

        $genotype="$ref/$alt[0]";
#    }
#elsif ($GT =~ m#^1/0#){
#        $zygosity = "het";
#        $alt_depth=(split(/,/,$AD))[0];
#        $genotype="$ref/$alt[0]";
    }elsif ($GT =~ m#^1/1#) {
        $zygosity =  "hom";
        
        if($AD=~/,/){
            $alt_depth=(split(/,/,$AD))[1];
        }else{
            $alt_depth=$AD;
        }

        $genotype="$alt[0]/$alt[0]";
     }elsif ($GT =~ m#^0/0#) {
        $zygosity = "ref";
        $genotype="$ref/$ref";
        $alt_depth="0";
    }elsif ($GT =~ m#^./.#) {
        $zygosity = "miss";
        $genotype = ".";
        $alt_depth="0";
    }elsif($GT=~m#2/1# || $GT=~m#1/2#){
        $zygosity="multalt";
        $genotype="$alt[1]/$alt[0]";
        my @AD=split(/,/,$AD);
        $alt_depth="$AD[2]/$AD[1]";
    }else{
        $zygosity="unkown";
        $genotype="NA";
        $alt_depth="NA";
    }
    

    return($zygosity,$genotype,$alt_depth);
}


sub SamZygJuge{
    ###based on cases and controls info to calculate mendel model,freq and pvalue 
    my ($sample_inf,$fam,$sample_fam,$ref,$alt,$format)=@_;
    ##@each sample genotype,$fam name, \@samples in name,$ref_allele,$mut_allele,vcf_format

   
    my @sample_info = @$sample_inf;
    my @sample_fam  = @$sample_fam;

    my (%hash_zyg,%gt,$out_alt);

### get each sample type 
    for(my $i=0;$i<=$#sample_fam;$i++){
        #next if not exists $group{$i};
        
        my($zygosity,$genotype,$alt_depth)=Zygosity($sample_info[$i],$ref,$alt);
     
        $out_alt.="$alt_depth;"; 
        $gt{$zygosity}=$genotype;
        
        if(not exists $hash_zyg{$zygosity}{$fam_n{$fam}{$sample_fam[$i]}}){
            $hash_zyg{$zygosity}{$fam_n{$fam}{$sample_fam[$i]}}=1;
        }else{
            $hash_zyg{$zygosity}{$fam_n{$fam}{$sample_fam[$i]}}++;
        }
    }
   


    my (@case_zyg,@control_zyg,$case_n,$control_n,@genotype);
    for my $type ('hom','het','ref','miss','multalt','unkown'){
        push @genotype,$gt{$type} if (exists $gt{$type});
       
        if(exists $hash_zyg{$type}{'2'}){
            push @case_zyg,$hash_zyg{$type}{'2'};
            $case_n+=$hash_zyg{$type}{'2'};
        }else{
            push @case_zyg,0;
            $case_n+=0;
        }

        if(exists $hash_zyg{$type}{'1'}){
            push @control_zyg,$hash_zyg{$type}{'1'};
            $control_n+=$hash_zyg{$type}{'1'};
        }else{
            push @control_zyg,0;
            $control_n+=0;
        }
    }
  
   my $nonemut_n=0; 
   for my $type ('ref','miss','unkown'){
        if(exists $hash_zyg{$type}{'2'}){
            $nonemut_n+=$hash_zyg{$type}{'2'};
        }
   }


   my($case_miss_rate,$control_miss_rate,$control_hom_rate,
      $control_het_rate,$model,$penetrance,$case_avgAltDp,$control_avgAltDp);
   
   $penetrance ||=0.8;

   $case_n==0?($case_miss_rate=0):($case_miss_rate=$case_zyg[4]/$case_n);

   $control_n>0 ? ($control_miss_rate=$control_zyg[4]/$control_n) : ($control_miss_rate=0);
   $control_n>0 ? ($control_hom_rate =$control_zyg[0]/$control_n) : ($control_hom_rate =0);
   $control_n>0 ? ($control_het_rate =$control_zyg[1]/$control_n) : ($control_het_rate =0);
   #print "$hash_zyg{'dp'}{'1'}\t$control_avgAltDp\t$control_miss_rate\t$control_hom_rate\t$control_het_rate\n";  
 
    
    if($case_zyg[0] >0 && $case_zyg[1]==0 && $case_zyg[5]==0 && $control_hom_rate<=1-$penetrance){
        #case only exists hom and miss,controls hom less than penetrance 
        if($case_zyg[0]==$case_n){#all cases are hom
            if($control_miss_rate==0){
                $model="AR:YY";
            }else{
                $model="AR:YN1";#control miss genotype,no mattter rate
            }
        }elsif($case_zyg[2]==0){#case no refgenotype,only have miss
            if($case_miss_rate <= 0.5){ #case miss genotype,allow 50% miss 
                if($control_miss_rate==0){
                    $model="AR:YN2";#case miss genotype
                }else{
                    $model="AR:YN12";#case and control miss genotype
                }
            }else{
                $model="AR:NN";
            }
        }else{
            $model="AR:NN";
        }
    }elsif($case_zyg[0]==0 && $control_zyg[0]==0 && $case_zyg[1]>0 && $case_zyg[5]==0 && $control_het_rate<=1-$penetrance){
        if($case_zyg[1]==$case_n){#all cases are het
            if($control_miss_rate==0){
                $model="AD:YY";
            }else{
                $model="AD:YN1";#control miss genotype,no mattter rate
            }
        }elsif($case_zyg[2]==0){#case no refgenotype,only have miss
            if($case_miss_rate <= 0.5){ #case miss genotype,allow 50% miss 
                if($control_miss_rate==0){
                    $model="AD:YN2";#case miss genotype
                }else{
                    $model="AD:YN12";#case and control miss genotype
                }
            }else{
                $model="AD:NN";
            }
        }else{
            $model="AD:NN";
        }
    }elsif($case_zyg[0]>0 && $case_zyg[1]>0){
        $model="conflic";
    }elsif($case_n==$nonemut_n){
        $model="NONEMUT";
    }else{
        $model="NONE";
    }
    

    #my $case_zyp_info = join (";",@case_zyg);#hom','het','ref','miss','multalt','unkown'
    #my $control_zyg_info = join (";",@control_zyg);#hom','het','ref','miss','multalt','unkown'
    my $case_zyp_info = join (";",@case_zyg[0..3]);#hom','het','ref','miss'
    my $control_zyg_info = join (";",@control_zyg[0..3]);#hom','het','ref','miss'
     
    my ($alles_info,$alt_freq,$twotailed_value) = calculate_freq_and_p(\@case_zyg,\@control_zyg);
    
    my $out_gt=join (",",@genotype);
    my $out_alt_dp;
        
    $model.=":".(join",",@genotype).":"."$case_n,$control_n:"; #."$case_zyp_info";
    return($model,$out_alt,$alles_info,$alt_freq,$twotailed_value);#mendel model;each sample alt depth;case_alt,case_ref:control_alt,case_ref,palue
}


sub getmMutSample{
  ## for  combine vcf,fetch mutated samples and its mutation type;
  my ($sample_inf,$fam,$sample_fam,$ref,$alt,$format)=@_;
  ##@each sample genotype,$fam name, \@samples in name,$ref_allele,$mut_allele,vcf_format
  
  my @sample_info = @$sample_inf;
  my @sample_fam  = @$sample_fam;

  my $fetch_sample="";
  my $fetch_sample_gt="";

  for(my $i=0;$i<=$#sample_fam;$i++){
    my($zygosity,$genotype,$alt_depth)=Zygosity($sample_info[$i],$ref,$alt);
    #zygosity is one of 'hom','het','ref','miss','multalt','unkown'
    foreach ('hom','het','multalt','unkown'){
      if ($zygosity eq $_){
        $fetch_sample.=$sample_fam[$i].",";
        $fetch_sample_gt.=$_.",";
        last;
      }
    }
  }
  return($fetch_sample,$fetch_sample_gt);
}


sub calculate_freq_and_p{
    #calculte  case and control freqency,fisher pvalue
    my ($case_zyg,$control_zyg)=@_;##case and control zyg array,#hom','het','ref','miss','multalt','unkown'
    my @case_zyg = @$case_zyg;
    my @control_zyg = @$control_zyg;

    my $case_alle_a = $case_zyg[0]*2 + $case_zyg[1]*1; #hom*2 + het*1
    my $case_alle_b = $case_zyg[2]*1; #ref*1
    my $control_alle_a = $control_zyg[0]*2 + $control_zyg[1]*1 ;#hom*2 + het*1
    my $control_alle_b = $control_zyg[2]*1 ;#ref*1 
   
    my ($case_freq,$control_freq);
    my $np1 = $case_alle_a + $case_alle_b;
    my $np2 = $control_alle_a + $control_alle_b;
    my $n1p = $case_alle_a + $control_alle_a;
    my $n2p = $case_alle_b + $control_alle_b; 

    ($np1 eq 0) ? ($case_freq = 0) : ($case_freq = sprintf("%.2f",$case_alle_a/$np1));
    ($np2 eq 0) ? ($control_freq = 0) : ($control_freq = sprintf("%.2f",$control_alle_a/$np2));
    
=cut###
          word2   ~word2
 word1    n11      n12 | n1p
~word1    n21      n22 | n2p
          --------------
          np1      np2   npp

$twotailed_value = calculateStatistic( n11=>$n11,
                                    n1p=>$n1p,
                                    np1=>$np1,
                                    npp=>$npp);
=cut###
    my $alles_info = "$case_alle_a,$case_alle_b:$control_alle_a,$control_alle_b";
    my $alt_freq = "$case_freq:$control_freq";
    my $twotailed_value = calculateStatistic( n11=>$case_alle_a,
                                              n1p=>$case_alle_a + $control_alle_a,
                                              np1=>$np1,
                                              npp=>$np1 + $np2
                                              );   
    #$twotailed_value = 0 if $twotailed_value eq ""; 
    return($alles_info,$alt_freq,$twotailed_value);
}

