#!/usr/bin/env ruby
require 'optparse'

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
     "sio: <http://semanticscience.org/resource/>",
     "up: <http://purl.uniprot.org/uniprot/>"
    ].each {|uri| print "@prefix #{uri} .\n"}
    print "\n"
  end
  module_function :prefixes

  Term2SO_gene = Hash[*[
  "IG_C_gene", "SO_0002123",
  "IG_C_pseudogene", "SO_0002100",
  "IG_D_gene", "SO_0002124",
  "IG_J_gene", "SO_0002125",
  "IG_J_pseudogene", "SO_0002101",
  "IG_pseudogene", "SO_0002098",
  "IG_V_gene", "SO_0002126",
  "IG_V_pseudogene", "SO_0002102",
  "lncRNA", "SO_0002127",
  "miRNA", "SO_0001265",
  "misc_RNA", "SO_0001263",
  "Mt_rRNA", "SO_0002363",
  "Mt_tRNA", "SO_0000088",
  "polymorphic_pseudogene", "SO_0001841",
  "processed_pseudogene", "SO_0000043",
  "protein_coding", "SO_0001217",
  "pseudogene", "SO_0000336",
  "ribozyme", "SO_0002181",
  "rRNA", "SO_0001637",
  "rRNA_pseudogene", "SO_0000777",
  "scaRNA", "SO_0002339",
  "scRNA", "SO_0001266",
  "snoRNA", "SO_0001267",
  "snRNA", "SO_0001268",
  "sRNA", "SO_0002342",
  "transcribed_processed_pseudogene", "SO_0002109",
  "transcribed_unitary_pseudogene", "SO_0002108",
  "transcribed_unprocessed_pseudogene", "SO_0002107",
  "translated_processed_pseudogene", "SO_0002105",
  "translated_unprocessed_pseudogene", "SO_0002106",
  "TR_C_gene", "SO_0002134",
  "TR_D_gene", "SO_0002135",
  "TR_J_gene", "SO_0002136",
  "TR_J_pseudogene", "SO_0002104",
  "TR_V_gene", "SO_0002137",
  "TR_V_pseudogene", "SO_0002103",
  "unitary_pseudogene", "SO_0001759",
  "unprocessed_pseudogene", "SO_0001760",
  "vault_RNA", "SO_0002358"]]

  Term2SO_transcript = Hash[*[
  "lncRNA", "SO_0001877",
  "miRNA", "SO_0000276",
  "misc_RNA", "SO_0000655",
  "Mt_rRNA", "SO_0002128",
  "Mt_tRNA", "SO_0002129",
  "nonsense_mediated_decay", "SO_0002114",
  "non_stop_decay", "SO_0002130",
  "processed_transcript", "SO_0001503",
  "protein_coding", "SO_0000234",
  "retained_intron", "SO_0000681",
  "ribozyme", "SO_0000374",
  "rRNA", "SO_0000252",
  "scaRNA", "SO_0002095",
  "scRNA", "SO_0000013",
  "snoRNA", "SO_0000275",
  "snRNA", "SO_0000274",
  "sRNA", "SO_0002247",
  "TEC", "SO_0002139",
  "vault_RNA", "SO_0000404"]]

  ENST = "http://rdf.ebi.ac.uk/resource/ensembl.transcript/"

  class TSV

    def initialize(file_name_1, file_name_2, file_name_3)
      @f_genes = open(file_name_1)
      @f_exons = open(file_name_2)
      @f_extl  = open(file_name_3)
      @gene_keys = set_keys(@f_genes)
      @exon_keys = set_keys(@f_exons)
      @extl_keys = set_keys(@f_extl)
      @gene_hash = {}
      @exon_hash = {}
      @extl_hash = {}
      @gene2transcripts = {}
      @transcript2protein = {}
      @transcript_hash = {}
      @gene2extls = {}
      @transcript2extls = {}
      STDERR.print "#{@extl_keys.join(' ')}\n"
      STDERR.print "parse external links ....\n"
      parse_external_links
      STDERR.print "parse ....\n"
      parse
      STDERR.print "rdf ....\n"
    end

    def set_keys(fh)
      ary = fh.gets.chomp.split("\t", -1)
      ary.map{|e| e.gsub(/[()]/, "").gsub(/[\/\-\s]/, "_").downcase.to_sym}
    end

    def parse_external_links()
      while line = @f_extl.gets
        vals = line.chomp.split("\t", -1)
        h = Hash[[@extl_keys, vals].transpose]
        @gene2extls[h[:gene_stable_id]] = h[:hgnc_id].sub("HGNC:", "") unless h[:hgnc_id] == ""
        unless h[:uniprotkb_swiss_prot_id] == ""
          if @transcript2extls.key?(h[:transcript_stable_id])
            @transcript2extls[h[:transcript_stable_id]] << h[:uniprotkb_swissprot_id]
          else
            @transcript2extls[h[:transcript_stable_id]] = [h[:uniprotkb_swissprot_id]]
          end
        end
        unless h[:uniprotkb_trembl_id] == ""
          if @transcript2extls.key?(h[:transcript_stable_id])
            @transcript2extls[h[:transcript_stable_id]] << h[:uniprotkb_trembl_id]
          else
            @transcript2extls[h[:transcript_stable_id]] = [h[:uniprotkb_trembl_id]]
          end
        end
#        STDERR.print "#{h}"
      end
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
                                            h[:strand].to_i]
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
                                                         h[:strand].to_i,
                                                         h[:transcript_name],
                                                         h[:transcript_type]
                                                       ]
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
        print "    a obo:#{Term2SO_gene[@gene_hash[gene_id][1]]} ;\n" if Term2SO_gene.key?(@gene_hash[gene_id][1])
        print "    rdfs:label \"#{@gene_hash[gene_id][0]}\" ;\n"
        print "    dcterms:identifier \"#{gene_id}\" ;\n"
        print "    dcterms:description \"#{@gene_hash[gene_id][2]}\" ;\n"
        print "    obo:RO_0002162 taxon:9606 ;\n"
        if @gene2extls.key?(gene_id)
          print "    rdfs:seeAlso <http://identifiers.org/ensembl/#{gene_id}> ,\n"
          print "                 hgnc:#{@gene2extls[gene_id]} ;\n"
        else
          print "    rdfs:seeAlso <http://identifiers.org/ensembl/#{gene_id}> ;\n"
        end
        print "    so:part_of <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
        print "    obo:BFO_0000050 <http://identifiers.org/hco/#{@gene_hash[gene_id][3]}#GRCh37> ;\n"
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
          print "enst:#{transcript_id} a term:#{@transcript_hash[transcript_id][6]} ;\n"
          print "    a obo:#{Term2SO_transcript[@transcript_hash[transcript_id][6]]} ;\n" if Term2SO_transcript.key?(@transcript_hash[transcript_id][6])
          print "    dcterms:identifier \"#{transcript_id}\" ;\n"
          print "    rdfs:label \"#{@transcript_hash[transcript_id][5]}\" ;\n"
          print "    obo:BFO_0000050 ens:#{gene_id} ;\n"
          print "    so:part_of ens:#{gene_id} ;\n"
          print "    so:transcribed_from ens:#{gene_id} ;\n"
          print "    so:translates_to ensp:#{@transcript2protein[transcript_id]} ;\n" unless @transcript2protein[transcript_id] == ""
          if @transcript2extls.key?(transcript_id)
            print "      rdfs:seeAlso #{@transcript2extls[transcript_id].map{|e| "up:#{e}"}.join(", ")} ;\n"
          end
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
            print "enst:#{transcript_id} obo:BFO_0000051 ense:#{exon[0]} .\n"
            print "enst:#{transcript_id} sio:SIO_000974 <#{ENST}#{transcript_id}#Exon_#{exon[3]}> .\n"
            print "<#{ENST}#{transcript_id}#Exon_#{exon[3]}> a sio:SIO_001261 ;\n"
            print "    sio:SIO_000628 ense:#{exon[0]} ;\n"
            print "    sio:SIO_000300 #{exon[3]} .\n"
            print "\n"
            unless exon_used.include?(exon[0])
            exon_used << exon[0]
            print "ense:#{exon[0]} a obo:SO_0000147 ;\n" # so:exon
            print "    rdfs:label \"#{exon[0]}\" ;\n"
            print "    dcterms:identifier \"#{exon[0]}\" ;\n"
            print "    so:part_of enst:#{transcript_id} ;\n"
            print "    obo:BFO_0000050 enst:#{transcript_id} ;\n"
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

def help
  print "Usage: ruby rdf_converter_ensembl_grch37.rb [options]\n"
  print "  -g, --gene path to the file for gene structures\n"
  print "  -e, --exon path to the file for exon structures\n"
  print "  -x, --xlink path to the file for cross links\n"
end

params = ARGV.getopts('g:e:d:o:x:h:', 'gene:', 'exon:', 'dir:', 'xlink:', 'grch:')
if (params["g"] || params["gene"]) && (params["e"] || params["exon"]) && (params["x"] || params["xlink"])
  if (params["g"])
    e = Ensembl::TSV.new(params["g"], params["e"], params["x"])
    Ensembl.prefixes
    e.rdf
    exit
  elsif (params["gene"])
    e = Ensembl::TSV.new(params["gene"], params["exon"], params["xlink"])
    Ensembl.prefixes
    e.rdf
    exit
  else
    help
    exit
  end
end



#fg = ARGV.shift
#fe = ARGV.shift
#e = Ensembl::TSV.new(fg, fe)
#Ensembl.prefixes
#e.rdf()

