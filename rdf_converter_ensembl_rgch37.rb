#!/usr/bin/env ruby

module Ensembl

  def prefixes
    ["m2r: <http://med2rdf.org/ontology/med2rdf#>",
     "cvo: <http://purl.jp/bio/10/clinvar/>",
     "faldo: <http://biohackathon.org/resource/faldo#>",
     "ens: <http://rdf.ebi.ac.uk/resource/ensembl/>",
     "ensp: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>",
     "enst: <http://rdf.ebi.ac.uk/resource/ensembl.transcript/>",
     "ense: <http://rdf.ebi.ac.uk/resource/ensembl.exon/>",
     "rdfs: <http://www.w3.org/2000/01/rdf-schema#>",
     "dcterms: <http://purl.org/dc/terms/>",
     "ensgido: <http://identifiers.org/ensembl/>",
     "term: <http://rdf.ebi.ac.uk/terms/ensembl/>",
     "identifiers: <http://identifiers.org/>",
     "obo: <http://purl.obolibrary.org/obo/>",
     "taxon: <http://identifiers.org/taxonomy/>",
     "hgnc: <http://identifiers.org/hgnc/>",
     "so: <http://purl.obolibrary.org/obo/so#>",
     "sio: <http://semanticscience.org/resource/>"
    ].each {|uri| print "@prefix #{uri} .\n"}
    print "\n"
  end
  module_function :prefixes

  Term2SO = Hash[*[
  "3prime_overlapping_ncrna", "",
  "antisense", "",
  "bidirectional_promoter_lncrna", "",
  "IG_C_gene", "",
  "IG_C_pseudogene", "",
  "IG_D_gene", "SO_0000510",
  "IG_J_gene", "",
  "IG_J_pseudogene", "",
  "IG_V_gene", "",
  "IG_V_pseudogene", "",
  "lincRNA", "SO_0001641",
  "LRG_gene", "",
  "macro_lncRNA", "",
  "miRNA", "SO_0001265",
  "misc_RNA", "SO_0000356",
  "Mt_rRNA", "",
  "Mt_tRNA", "SO_0000088",
  "non_coding", "",
  "nonsense_mediated_decay", "SO_0001621",
  "non_stop_decay", "",
  "polymorphic_pseudogene", "SO_0000336",
  "processed_pseudogene", "",
  "processed_transcript", "SO_0001503",
  "protein", "",
  "protein_coding", "SO_0001217",
  "pseudogene", "SO_0000336",
  "retained_intron", "SO_0000681",
  "ribozyme", "",
  "rRNA", "SO_0001637",
  "scaRNA", "",
  "sense_intronic", "",
  "sense_overlapping", "",
  "snoRNA", "SO_0001267",
  "snRNA", "SO_0001268",
  "sRNA", "",
  "TEC", "",
  "transcribed_processed_pseudogene", "",
  "transcribed_unitary_pseudogene", "",
  "transcribed_unprocessed_pseudogene", "",
  "translated_unprocessed_pseudogene", "",
  "TR_C_gene", "SO_0000478",
  "TR_D_gene", "",
  "TR_J_gene", "SO_0000470",
  "TR_J_pseudogene", "",
  "TR_V_gene", "SO_0000466",
  "TR_V_pseudogene", "",
  "unitary_pseudogene", "",
  "unprocessed_pseudogene", "",
  "vaultRNA", "" ]]

  ENST = "http://rdf.ebi.ac.uk/resource/ensembl.transcript/"

  class TSV

    def initialize(file_name_1, file_name_2)
      @f_genes = open(file_name_1)
      @f_exons = open(file_name_2)
      @gene_keys = set_keys(@f_genes)
      @exon_keys = set_keys(@f_exons)
      @gene_hash = {}
      @exon_hash = {}
      @gene2transcripts = {}
      @transcript2protein = {}
      @transcript_hash = {}
      parse
    end

    def set_keys(fh)
      ary = fh.gets.chomp.split("\t", -1)
      ary.map{|e| e.gsub(/[()]/, "").gsub(/[\/\s]/, "_").downcase.to_sym}
    end

    def parse()
      while line = @f_genes.gets
        vals = line.chomp.split("\t", -1)
        h = Hash[[@gene_keys, vals].transpose]
        unless @gene_hash.key?(h[:gene_stable_id])
          @gene_hash[h[:gene_stable_id]] = [h[:gene_name],
                                            h[:gene_type],
                                            h[:gene_description],
                                            h[:chromosome_scaffold_name],
                                            h[:gene_start_bp].to_i,
                                            h[:gene_end_bp].to_i,
                                            h[:strand].to_i,
                                            h[:hgnc_id]]
        end
        if @gene2transcripts.key?(h[:gene_stable_id])
          unless @gene2transcripts[h[:gene_stable_id]].include?(h[:transcript_stable_id])
            @gene2transcripts[h[:gene_stable_id]] << h[:transcript_stable_id]
          end
        else
          @gene2transcripts[h[:gene_stable_id]] = [h[:transcript_stable_id]]
        end
        unless @transcript2protein.key?(h[:transcript_stabel_id])
          @transcript2protein[h[:transcript_stable_id]] = h[:protein_stable_id]
        end
        unless @transcript_hash.key?(h[:transcript_stable_id])
          @transcript_hash[h[:transcript_stable_id]] = [[h[:exon_stable_id]],
                                                         h[:chromosome_scaffold_name],
                                                         h[:transcript_start_bp].to_i,
                                                         h[:transcript_end_bp].to_i,
                                                         h[:strand].to_i]
        else
          @transcript_hash[h[:transcript_stable_id]][0] << h[:exon_stable_id]
        end
      end

      while line = @f_exons.gets
        vals = line.chomp.split("\t", -1)
        h = Hash[[@exon_keys, vals].transpose]
        unless @exon_hash.key?(h[:transcript_stable_id])
          @exon_hash[h[:transcript_stable_id]] = [[h[:exon_stable_id],
                                                   h[:exon_region_start_bp].to_i,
                                                   h[:exon_region_end_bp].to_i,
                                                   h[:exon_rank_in_transcript].to_i]]
        else
          @exon_hash[h[:transcript_stable_id]] << [h[:exon_stable_id],
                                                   h[:exon_region_start_bp].to_i,
                                                   h[:exon_region_end_bp].to_i,
                                                   h[:exon_rank_in_transcript].to_i]
        end
      end
    end

    def rdf()
      exon_used = []
      @gene_hash.keys.each do |gene_id|
        print "ens:#{gene_id} a term:#{@gene_hash[gene_id][1]} ;\n"
        print "    a obo:#{Term2SO[@gene_hash[gene_id][1]]} ;\n" unless Term2SO[@gene_hash[gene_id][1]] == ""
        print "    rdfs:label \"#{@gene_hash[gene_id][0]}\" ;\n"
        print "    dcterms:identifier \"#{gene_id}\" ;\n"
        print "    obo:RO_0002162 taxon:9606 ;\n"
        print "    rdfs:seeAlso <http://identifiers.org/ensembl/#{gene_id}> ,\n"
        print "                 hgnc:HGNC_#{@gene_hash[gene_id][7]} ;\n"
        print "    so:part_of <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
        print "    faldo:location [\n"
        print "        a faldo:Region ;\n"
        print "        faldo:begin [\n"
        print "            a faldo:ExactPosition ;\n"
        if @gene_hash[gene_id][6] == 1
          print "            a faldo:ForwardStrandPosition ;\n"
        else
          print "            a faldo:ReverseStrandPosition ;\n"
        end
        print "            faldo:position #{@gene_hash[gene_id][4]} ;\n"
        #  print "            faldo:reference hco:#{gene_hash[gene_id][3]}\\#GRCh37\n"
        print "            faldo:reference <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
        print "        ] ;\n"
        print "        faldo:end [\n"
        print "            a faldo:ExactPosition ;\n"
        if @gene_hash[gene_id][6] == 1
          print "            a faldo:ForwardStrandPosition ;\n"
        else
          print "            a faldo:ReverseStrandPosition ;\n"
        end
        print "            faldo:position #{@gene_hash[gene_id][5]} ;\n"
        #  print "            faldo:reference hco:#{gene_hash[gene_id][3]}\\#GRCh37\n"
        print "            faldo:reference <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
        print "        ]\n"
        print "    ] .\n"
        print "\n"
        @gene2transcripts[gene_id].each do |transcript_id|
          print "enst:#{transcript_id} a term:#{@gene_hash[gene_id][1]} ;\n"
          print "    a ens:#{Term2SO[@gene_hash[gene_id][1]]} ;\n" unless Term2SO[@gene_hash[gene_id][1]] == ""
          print "    dcterms:identifier \"#{transcript_id}\" ;\n"
          print "    so:part_of ens:#{gene_id} ;\n"
          print "    so:transcribed_from ens:#{gene_id} ;\n"
          print "    so:translates_to ensp:#{@transcript2protein[transcript_id]} ;\n" unless @transcript2protein[transcript_id] == ""
          print "    faldo:location [\n"
          print "        a faldo:Region ;\n"
          print "        faldo:begin [\n"
          print "            a faldo:ExactPosition ;\n"
          if @gene_hash[gene_id][6] == 1
            print "            a faldo:ForwardStrandPosition ;\n"
          else
            print "            a faldo:ReverseStrandPosition ;\n"
          end
          print "            faldo:position #{@transcript_hash[transcript_id][2]} ;\n"
          print "            faldo:reference <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
          print "        ] ;\n"
          print "        faldo:end [\n"
          print "            a faldo:ExactPosition ;\n"
          if @gene_hash[gene_id][6] == 1
            print "            a faldo:ForwardStrandPosition ;\n"
          else
            print "            a faldo:ReverseStrandPosition ;\n"
          end
          print "            faldo:position #{@transcript_hash[transcript_id][3]} ;\n"
          print "            faldo:reference <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
          print "        ]\n"
          print "    ] .\n"
          print "\n"
          @exon_hash[transcript_id].each do |exon|
            print "enst:#{transcript_id} so:has_part ense:#{exon[0]} .\n"
            print "enst:#{transcript_id} sio:SIO_000974 <#{ENST}#{transcript_id}#Exon_#{exon[3]}> .\n"
            print "<#{ENST}#{transcript_id}#Exon_#{exon[3]}> a sio:SIO_001261 ;\n"
            print "    sio:SIO_000628 ense:#{exon[0]} ;\n"
            print "    sio:SIO_000300 #{exon[3]} .\n"
            print "\n"
            unless exon_used.include?(exon[0])
            exon_used << exon[0]
            print "ense:#{exon[0]} a obo:SO_0000147 ;\n" # so:exon
            print "    rdfs:label \"#{exon[0]}\" ;\n"
            print "    dcterms:identifiers \"#{exon[0]}\" ;\n"
            print "    so:part_of enst:#{transcript_id} ;\n"
            print "    faldo:location [\n"
            print "        a faldo:Region ;\n"
            print "        faldo:begin [\n"
            print "            a faldo:ExactPosition ;\n"
            if @gene_hash[gene_id][6] == 1
              print "            a faldo:ForwardStrandPosition ;\n"
            else
              print "            a faldo:ReverseStrandPosition ;\n"
            end
            print "            faldo:position #{exon[1]} ;\n"
            print "            faldo:reference <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
            print "        ] ;\n"
            print "        faldo:end [\n"
            print "            a faldo:ExactPosition ;\n"
            if @gene_hash[gene_id][6] == 1
              print "            a faldo:ForwardStrandPosition ;\n"
            else
              print "            a faldo:ReverseStrandPosition ;\n"
            end
            print "            faldo:position #{exon[2]} ;\n"
            print "            faldo:reference <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
            print "        ] \n"
            print "    ] .\n"
            print "\n"
            end
          end
        end
        print "<http://identifiers.org/ensembl/#{gene_id}> a identifiers:ensembl .\n"
        print "\n"
      end
    end # end of rdf()
  end
end # end of Module

fg = ARGV.shift
fe = ARGV.shift
e = Ensembl::TSV.new(fg, fe)
Ensembl.prefixes
e.rdf()

=begin

term2so = Hash[*[
  "3prime_overlapping_ncrna", "",
  "antisense", "",
  "bidirectional_promoter_lncrna", "",
  "IG_C_gene", "",
  "IG_C_pseudogene", "",
  "IG_D_gene", "SO_0000510",
  "IG_J_gene", "",
  "IG_J_pseudogene", "",
  "IG_V_gene", "",
  "IG_V_pseudogene", "",
  "lincRNA", "SO_0001641",
  "LRG_gene", "",
  "macro_lncRNA", "",
  "miRNA", "SO_0001265",
  "misc_RNA", "SO_0000356",
  "Mt_rRNA", "",
  "Mt_tRNA", "SO_0000088",
  "non_coding", "",
  "nonsense_mediated_decay", "SO_0001621",
  "non_stop_decay", "",
  "polymorphic_pseudogene", "SO_0000336",
  "processed_pseudogene", "",
  "processed_transcript", "SO_0001503",
  "protein", "",
  "protein_coding", "SO_0001217",
  "pseudogene", "SO_0000336",
  "retained_intron", "SO_0000681",
  "ribozyme", "",
  "rRNA", "SO_0001637",
  "scaRNA", "",
  "sense_intronic", "",
  "sense_overlapping", "",
  "snoRNA", "SO_0001267",
  "snRNA", "SO_0001268",
  "sRNA", "",
  "TEC", "",
  "transcribed_processed_pseudogene", "",
  "transcribed_unitary_pseudogene", "",
  "transcribed_unprocessed_pseudogene", "",
  "translated_unprocessed_pseudogene", "",
  "TR_C_gene", "SO_0000478",
  "TR_D_gene", "",
  "TR_J_gene", "SO_0000470",
  "TR_J_pseudogene", "",
  "TR_V_gene", "SO_0000466",
  "TR_V_pseudogene", "",
  "unitary_pseudogene", "",
  "unprocessed_pseudogene", "",
  "vaultRNA", "" ]]


print "@prefix m2r: <http://med2rdf.org/ontology/med2rdf#> .\n"
print "@prefix cvo: <http://purl.jp/bio/10/clinvar/> .\n"
print "@prefix faldo: <http://biohackathon.org/resource/faldo#> .\n"
print "@prefix ens: <http://rdf.ebi.ac.uk/resource/ensembl/> .\n"
print "@prefix ensp: <http://rdf.ebi.ac.uk/resource/ensembl.protein/> .\n"
print "@prefix enst: <http://rdf.ebi.ac.uk/resource/ensembl.transcript/> .\n"
print "@prefix ense: <http://rdf.ebi.ac.uk/resource/ensembl.exon/> .\n"
print "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .\n"
print "@prefix dcterms: <http://purl.org/dc/terms/> .\n"
print "@prefix ensgido: <http://identifiers.org/ensembl/> .\n"
print "@prefix term: <http://rdf.ebi.ac.uk/terms/ensembl/> .\n"
print "@prefix identifiers: <http://identifiers.org/> .\n"
print "@prefix obo: <http://purl.obolibrary.org/obo/> .\n"
print "@prefix so: <http://purl.obolibrary.org/obo/so#> .\n"
print "@prefix sio: <http://semanticscience.org/resource/> .\n"
print "\n"

ENST = "http://rdf.ebi.ac.uk/resource/ensembl.transcript/"

f_genes = open(ARGV.shift) # genes
f_exons = open(ARGV.shift) # structures
keys = f_genes.gets.chomp.split("\t", -1).map{|e| e.gsub("(", "").gsub(")","").gsub("/", " ").gsub(" ", "_").downcase}
exon_keys = f_exons.gets.chomp.split("\t", -1).map{|e| e.gsub("(", "").gsub(")","").gsub("/", " ").gsub(" ", "_").downcase}

gene_hash = {}
gene2transcripts = {}
transcript_hash = {}

while line = f_genes.gets
  vals = line.chomp.split("\t", -1)
  h = Hash[[keys, vals].transpose]
  unless gene_hash.key?(h["gene_stable_id"])
    gene_hash[h["gene_stable_id"]] = [h["gene_name"],
                                      h["gene_type"],
                                      h["gene_description"],
                                      h["chromosome_scaffold_name"],
                                      h["gene_start_bp"].to_i,
                                      h["gene_end_bp"].to_i,
                                      h["strand"].to_i,
                                      h["protein_stable_id"],
                                      h["hgnc_id"]]
  end
  unless gene2transcripts.key?(h["gene_stable_id"])
    gene2transcripts[h["gene_stable_id"]] = [h["transcript_stable_id"]]
  else
    gene2transcripts[h["gene_stable_id"]] << h["transcript_stable_id"]
  end
  unless transcript_hash.key?(h["transcript_stable_id"])
    transcript_hash[h["transcript_stable_id"]] = [[h["exon_stable_id"]],
                                                  h["chromosome_scaffold_name"],
                                                  h["transcript_start_bp"].to_i,
                                                  h["transcript_end_bp"].to_i,
                                                  h["strand"].to_i]
  else
    transcript_hash[h["transcript_stable_id"]][0] << h["exon_stable_id"]
  end
end

exon_hash = {}

while line = f_exons.gets
  vals = line.chomp.split("\t", -1)
  h = Hash[[exon_keys, vals].transpose]
  unless exon_hash.key?(h["transcript_stable_id"])
    exon_hash[h["transcript_stable_id"]] = [[h["exon_stable_id"],
                                             h["exon_region_start_bp"].to_i,
                                             h["exon_region_end_bp"].to_i,
                                             h["exon_rank_in_transcript"].to_i]]
  else
    exon_hash[h["transcript_stable_id"]] << [h["exon_stable_id"],
                                             h["exon_region_start_bp"].to_i,
                                             h["exon_region_end_bp"].to_i,
                                             h["exon_rank_in_transcript"].to_i]

  end
end


gene_hash.keys.each do |gene_id|
  print "ens:#{gene_id} a term:#{gene_hash[gene_id][1]} ;\n"
  print "    a obo:#{term2so[gene_hash[gene_id][1]]} ;\n" unless term2so[gene_hash[gene_id][1]] == ""
  print "    rdfs:label \"#{gene_hash[gene_id][0]}\" ;\n"
  print "    dcterms:identifier \"#{gene_id}\" ;\n"
  print "    rdfs:seeAlso <http://identifiers.org/ensembl/#{gene_id}> ;\n"
  print "    so:part_of <http://identifiers.org/hco/#{gene_hash[gene_id][3]}#GRCh37> ;\n"
  print "    faldo:location [\n"
  print "        a faldo:Region ;\n"
  print "        faldo:begin [\n"
  print "            a faldo:ExactPosition ;\n"
  if gene_hash[gene_id][6] == 1
  print "            a faldo:ForwardStrandPosition ;\n"
  else
  print "            a faldo:ReverseStrandPosition ;\n"
  end
  print "            faldo:position #{gene_hash[gene_id][4]} ;\n"
#  print "            faldo:reference hco:#{gene_hash[gene_id][3]}\\#GRCh37\n"
  print "            faldo:reference <http://identifiers.org/hco/#{gene_hash[gene_id][3]}#GRCh37> ;\n"
  print "        ] ;\n"
  print "        faldo:end [\n"
  print "            a faldo:ExactPosition ;\n"
  if gene_hash[gene_id][6] == 1
  print "            a faldo:ForwardStrandPosition ;\n"
  else
  print "            a faldo:ReverseStrandPosition ;\n"
  end
  print "            faldo:position #{gene_hash[gene_id][5]} ;\n"
#  print "            faldo:reference hco:#{gene_hash[gene_id][3]}\\#GRCh37\n"
  print "            faldo:reference <http://identifiers.org/hco/#{gene_hash[gene_id][3]}#GRCh37> ;\n"
  print "        ]\n"
  print "    ] .\n"
  print "\n"
  gene2transcripts[gene_id].each do |transcript_id|
    print "enst:#{transcript_id} a term:#{gene_hash[gene_id][1]} ;\n"
    print "    a ens:#{term2so[gene_hash[gene_id][1]]} ;\n" unless term2so[gene_hash[gene_id][1]] == ""
    print "    dcterms:identifier \"#{transcript_id}\" ;\n"
    print "    so:part_of ens:#{gene_id} ;\n"
    print "    so:transcribed_from ens:#{gene_id} ;\n"
    print "    so:translates_to ensp:#{transcript_hash[transcript_id][6]} ;\n" unless transcript_hash[transcript_id][6] == ""
    print "    faldo:location [\n"
    print "        a faldo:Region ;\n"
    print "        faldo:begin [\n"
    print "            a faldo:ExactPosition ;\n"
    if gene_hash[gene_id][6] == 1
    print "            a faldo:ForwardStrandPosition ;\n"
    else
    print "            a faldo:ReverseStrandPosition ;\n"
    end
    print "            faldo:position #{transcript_hash[transcript_id][2]} ;\n"
    print "            faldo:reference <http://identifiers.org/hco/#{gene_hash[gene_id][3]}#GRCh37> ;\n"
    print "        ] ;\n"
    print "        faldo:end [\n"
    print "            a faldo:ExactPosition ;\n"
    if gene_hash[gene_id][6] == 1
    print "            a faldo:ForwardStrandPosition ;\n"
    else
    print "            a faldo:ReverseStrandPosition ;\n"
    end
    print "            faldo:position #{transcript_hash[transcript_id][3]} ;\n"
    print "            faldo:reference <http://identifiers.org/hco/#{gene_hash[gene_id][3]}#GRCh37> ;\n"
    print "        ]\n"
    print "    ] .\n"
    print "\n"
    exon_hash[transcript_id].each do |exon|
      print "enst:#{transcript_id} so:has_part ense:#{exon[0]} .\n"
      print "enst:#{transcript_id} sio:SIO_000974 <#{ENST}#{transcript_id}#Exon_#{exon[3]}> .\n"
      print "<#{ENST}#{transcript_id}#Exon_#{exon[3]}> a sio:SIO_001261 ;\n"
      print "    sio:SIO_000628 ense:#{exon[0]} ;\n"
      print "    sio:SIO_000300 #{exon[3]} .\n"
      print "ense:#{exon[0]} a obo:SO_0000147 ;\n" # so:exon
      print "    rdfs:label \"#{exon[0]}\" ;\n"
      print "    dcterms:identifiers \"#{exon[0]}\" ;\n"
      print "    so:part_of enst:#{transcript_id} ;\n"
      print "    faldo:location [\n"
      print "        a faldo:Region ;\n"
      print "        faldo:begin [\n"
      print "            a faldo:ExactPosition ;\n"
      if gene_hash[gene_id][6] == 1
      print "            a faldo:ForwardStrandPosition ;\n"
      else
      print "            a faldo:ReverseStrandPosition ;\n"
      end
      print "            faldo:position #{exon[1]} ;\n"
      print "            faldo:reference <http://identifiers.org/hco/#{gene_hash[gene_id][3]}#GRCh37> ;\n"
      print "        ] ;\n"
      print "        faldo:end [\n"
      print "            a faldo:ExactPosition ;\n"
      if gene_hash[gene_id][6] == 1
      print "            a faldo:ForwardStrandPosition ;\n"
      else
      print "            a faldo:ReverseStrandPosition ;\n"
      end
      print "            faldo:position #{exon[2]} ;\n"
      print "            faldo:reference <http://identifiers.org/hco/#{gene_hash[gene_id][3]}#GRCh37> ;\n"
      print "        ] \n"
      print "    ] .\n"
      print "\n"
    end
  end
  print "<http://identifiers.org/ensembl/#{gene_id}> a identifiers:ensembl .\n"
  print "\n"
end

=end
