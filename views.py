...


class SearchVariant(LoginRequiredMixin, View):
    template_name = 'lab/search_variant.html'

    def __init__ (self):
        
        # var init
        self.genome_version = ""
        self.filtred_variants =  []
        self.filtred_report = []

        # var init outsource_variants
        self.filtred_outsource_variant = [] 
        self.filtred_outsource_reports = []
        self.outsource_variant_atributes_dict = {}


    def genome_name_convert(self, chrom, position, ref, alt):
        if len(ref) == 1 and len(alt) == 1 :   
            genome_name = (f"chr{chrom}:g.{position}{ref}>{alt}")

        elif len(ref)>len(alt) and len(alt) == 1:
            position = int(position)+1
            del_size = len(ref)-2 + position

            if len(ref) == 2:
                genome_name = (f"chr{chrom}:g.{position}del")
            else:
                genome_name = (f"chr{chrom}:g.{position}_{del_size}del")   

        elif len(ref)<len(alt) and len(ref)== 1:
            ins_size = int(position)+1
            genome_name = (f"chr{chrom}:g.{position}_{ins_size}ins{alt[1:]}")

        elif len(ref)>=len(alt) or len(ref)<len(alt):
            position = int(position)+1
            if len(ref)>=len(alt):
                ins_del_size = len(ref)-2 + position
            elif len(ref)<len(alt):
                ins_del_size = len(alt)-2 + position 

            genome_name = (f"chr{chrom}:g.{position}_{ins_del_size}delins{alt[1:]}")

        return genome_name


    #https://grch37.rest.ensembl.org/variant_recoder/human/{rsid_or_variant_c}?fields=hgvsg&vcf_string=1&content-type=application/json
    def rsid_or_c_to_g(self, rsid_or_variant_c, genome_version):
        if genome_version == "GRCh38":
            server = "https://rest.ensembl.org" 
        elif genome_version == "GRCh37": 
            server = "https://grch37.rest.ensembl.org"

        ext = f"/variant_recoder/human/{rsid_or_variant_c}?"
        params ="fields=hgvsg&vcf_string=1&content-type=application/json"

        response = requests.get(server+ext+params)
        #response.raise_for_status()  # Raise an exception for non-OK responses
        
        if not response.ok:
            return []

        response_json = response.json()
        vcf_strings = []

        for key in response_json[0]:
            try:
                if 'vcf_string' in response_json[0][key]:
                    vcf_strings.append(response_json[0][key]['vcf_string'][0])
            
            except KeyError:
                pass
                
        return vcf_strings   

    def get_rsid_transcript_pos_ref_alt_gene(self,searchquery):
        rsid, transcript, gene, pos_ref_alt = "", "", "", ""
        
        #regex to get all digits after rs or all digits after NM_ or all digits and then everything after c.
        regex = r"(rs\d+)|(NM_\d+)|(c\.\d+.*)"
        matches = re.findall(regex, searchquery)
        for match in matches:
            if match[0]:
                rsid = match[0]
            elif match[1]:
                transcript = match[1]
            elif match[2]:
                pos_ref_alt = match[2]
            if not transcript and pos_ref_alt:
                gene = searchquery.replace(pos_ref_alt,"")
        return rsid, transcript, pos_ref_alt, gene


    def get_searchquery_converted_list(self, rsid, transcript, pos_ref_alt, gene, genome_version):
        if rsid:
            # convert rsid to a list of vcf_strings
            searchquery_converted_list = self.rsid_or_c_to_g(rsid, genome_version)
            return searchquery_converted_list

        if transcript:
            # convert NM_ :pos_ref_alt to list of vcf_strings
            searchquery_converted_list = self.rsid_or_c_to_g(transcript + ":" + pos_ref_alt, genome_version) 
            return searchquery_converted_list
        
        if gene:
            filtered_transcripts = Transcript.objects.filter(gene__hgnc_symbol=gene).values_list('refseq_mrna_noversion', flat=True)

            # put all variants in format NM_:pos_alt_ref
            searchquery_converted_list = [filtred_transcript + ":" + pos_ref_alt for filtred_transcript in filtered_transcripts]

            # convert all NM_variants to a list of vcf_strings
            searchquery_converted_list = [self.rsid_or_c_to_g(NM_variant, genome_version) for NM_variant in searchquery_converted_list]

            # if the rsid_or_c_to_g(NM_variant) return a nested_list we need to convert every element to one list
            searchquery_converted_list = [element for sublist in searchquery_converted_list for element in sublist]
            return searchquery_converted_list

    
    def get_genomic_coord_list_or_chrom_position_list(self, searchquery_converted_list):
        # chrom-pos-ref-alt
        genomic_coord_list = []

        # chrom-pos
        chrom_position_list = []

        # check all variants in searchquery_converted_list
        for variant_vcf in searchquery_converted_list:
            variant_atributes_list = variant_vcf.split("-")

            if len(variant_atributes_list) == 4:
                chrom, pos, ref, alt = variant_atributes_list[0],variant_atributes_list[1],variant_atributes_list[2],variant_atributes_list[3]
                genomic_coord = self.genome_name_convert(chrom, pos, ref, alt)
                genomic_coord_list.append( genomic_coord)

            elif len(variant_atributes_list) == 2:
                chrom_position_list.append(variant_vcf)
        
        return genomic_coord_list, chrom_position_list         
        

    
    def get_filtred_outsource_variant(self, rsid, transcript, gene, pos_ref_alt):        
        # it returns [] when the qsearchquery is in (chrom-pos-ref-alt) format othewise ( if i dont have this ini var) will return None
        filtred_outsource_variant = []
        
        if rsid:
            filtred_outsource_variant = RunSampleVariantOutsourcedData.objects.filter(dbsnp_id=rsid)
            return filtred_outsource_variant

        if transcript:
            filtred_outsource_variant = RunSampleVariantOutsourcedData.objects.filter(datasource_transcript__transcript__refseq_mrna_noversion=transcript, hgvs_nomenclature__icontains=pos_ref_alt)
            return filtred_outsource_variant
        
        if gene:
            filtred_outsource_variant = RunSampleVariantOutsourcedData.objects.filter(gene__hgnc_symbol=gene, hgvs_nomenclature__icontains=pos_ref_alt)
            return filtred_outsource_variant
        return filtred_outsource_variant

    def get_outsource_variant_atributes_dict(self, filtred_outsource_variant):
        # mesmo que o resultado do filtred_outsource_variant sejam varias variantes  
        # a variante é só uma, pois, a query é do tipo NM_:c. ou gene:c ou rs_id

        outsource_variant_atributes_dict = {}
        filtred_outsource_reports = []
        sample_names_list = []

        for outsource_variant in filtred_outsource_variant:
            chromossome_position = outsource_variant.position
            regex = r"chr(\d{2}|\d|X|Y):(\d+.*)"
            matches = re.findall(regex, chromossome_position, re.IGNORECASE)
            for match in matches:
                chromosome = match[0] if match[0] else "-"
                position = match[1] if match[1] else "-"

                if 'chromosome' not in outsource_variant_atributes_dict or not outsource_variant_atributes_dict['chromosome']:
                    outsource_variant_atributes_dict['chromosome'] = chromosome
                if 'position' not in outsource_variant_atributes_dict or not outsource_variant_atributes_dict['position']:
                    outsource_variant_atributes_dict['position'] = position

            gene = outsource_variant.gene.hgnc_symbol if outsource_variant.gene else "-"
            dbsnp = outsource_variant.dbsnp_id if outsource_variant.dbsnp_id else "-"
            hgvsc = outsource_variant.hgvs_nomenclature if outsource_variant.hgvs_nomenclature else "-"
            sample_names_list.append(outsource_variant.sample.sample_name)

            if not outsource_variant_atributes_dict.get('gene'):
                outsource_variant_atributes_dict['gene'] = gene
            if not outsource_variant_atributes_dict.get('dbsnp'):
                outsource_variant_atributes_dict['dbsnp'] = dbsnp
            if not outsource_variant_atributes_dict.get('hgvsc'):
                outsource_variant_atributes_dict['hgvsc'] = hgvsc
                

        outsource_variant_atributes_dict['samples'] = sample_names_list
        filtred_outsource_reports = Report.objects.filter(sample__sample_name__in=sample_names_list)

        return outsource_variant_atributes_dict, filtred_outsource_reports


    def get(self, request):
        searchquery = request.GET.get('searchvariant')
        
        if not searchquery:
            return render(request, self.template_name)
        
        regex = r"(GRCh3[7-8])"
        match = re.search(regex, searchquery)
        if not match:
            return render(request, self.template_name)
        
        self.genome_version = match.group(0)

        searchquery = searchquery.replace(self.genome_version, '') 
        


        # remove all spaces in searchquery
        searchquery = searchquery.replace(' ', '') 


        # remove all ":" in searchquery
        searchquery = searchquery.replace(':','')


        # get rsid transcript pos_ref_alt and gene from searchquery input
        rsid, transcript, pos_ref_alt, gene = self.get_rsid_transcript_pos_ref_alt_gene(searchquery)
        
        
        # if its in format chrom-pos-ref-alt or chrom-pos or error
        if not rsid and not transcript and not pos_ref_alt and not gene:
            searchquery_converted_list = [searchquery]
        

        # all formats that contains rsid, transcript, pos_ref_alt, gene
        else: searchquery_converted_list =  self.get_searchquery_converted_list(rsid, transcript, pos_ref_alt, gene, self.genome_version)


        # get a list of all genomic_coord_list if the searchquery_converted_list is in formart [chrom-pos-ref-alt...]
        # or get a list of all chrom_position_list is the searchquery_converted_list is in format [chrom-pos...]
        # if the searchquery_converted_list is in another format return [] per list
        genomic_coord_list, chrom_position_list  = self.get_genomic_coord_list_or_chrom_position_list(searchquery_converted_list)


        # get genomic filtred_variants (from genetic_variants table)
        # if dosent exist get the filtred_outsource_variant (from outsource_variants table) and then the outsource_variant_atributes_dict
        if genomic_coord_list:
            print(genomic_coord_list)
            self.filtred_variants = Variant.objects.filter(id__in=genomic_coord_list)
            if not self.filtred_variants:
                self.filtred_outsource_variant = self.get_filtred_outsource_variant(rsid, transcript, gene, pos_ref_alt)
                self.outsource_variant_atributes_dict, self.filtred_outsource_reports = self.get_outsource_variant_atributes_dict(self.filtred_outsource_variant)
            else: self.filtred_report = VariantClassification.objects.filter(variant__in= self.filtred_variants)


        # get only chrom_position filtred_variants (from genetic_variants table) 
        # if dosent exist range the query position to +10-10 bp  
        for variant_chr_pos in chrom_position_list:
            variant_atributes_list = variant_chr_pos.split("-")
            self.filtred_variants = Variant.objects.filter(chrom=variant_atributes_list[0],position=variant_atributes_list[1])

            if not self.filtred_variants:
                pos = variant_atributes_list[1]

                #Vai até +10
                position_range = range(int(pos)-10, int(pos)+11) 
                self.filtred_variants = Variant.objects.filter(chrom=variant_atributes_list[0], position__in=position_range) 
                
                #Position ajust alert
                messages.warning(request,"Variant not found in Ophidx Database!")

            self.filtred_report = VariantClassification.objects.filter(variant__in=self.filtred_variants)
        

        #Catch errors
        if self.filtred_outsource_variant:
            messages.success(request,"Variant found in Outsource Database!")
        
        # Variant not found in VEP! 
        if not searchquery_converted_list: 
            messages.error(request,"Variant not found in VEP!")

        # Bad terms   
        elif self.filtred_variants == []:
            messages.info(request, "Bad Terms!")
            
        
        context = { 
            'variants': self.filtred_variants,
            'reports': self.filtred_report,
            'searchquery': searchquery,
            'searchquery_converted':searchquery_converted_list,

            'outsource_variants': self.filtred_outsource_variant,
            'outsource_reports': self.filtred_outsource_reports,
            'outsource_variant_atributes' : self.outsource_variant_atributes_dict
        }

        return render(request, self.template_name, context)

