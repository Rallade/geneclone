package clone.gene;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class DNA{
	
	private int GENEPERCHROM = 3; //3 genes per chromosome
	private int CHROMSIZE = 50;
	private int ploidy;
	
	public DNA(int chromosomes, int homologuous){
		Chromosome c;
		for(int i = 0; i < chromosomes; i++){
			c = new Chromosome(homologuous);
			c.setRand(CHROMSIZE, GENEPERCHROM);
			dna.add(c);
		}
		ploidy = homologuous;
	}
	
	public DNA(ArrayList<ArrayList<ArrayList<Character>>> sequence) { //possiblity to include polysomy
		Nucleotide n;
		ArrayList<Nucleotide> chromatidN = new ArrayList<>();
		ArrayList<ArrayList<Nucleotide>> chromosomeN = new ArrayList<>();
		Chromosome chrom;
		
		for(ArrayList<ArrayList<Character>> chromosomeC : sequence){
			for(ArrayList<Character> chromatidC : chromosomeC){
				for(char c : chromatidC){
					n = new Nucleotide(c);
					chromatidN.add(n);
				}
				chromosomeN.add(chromatidN);
				chromatidN = new ArrayList<>();
			}
			chrom = new Chromosome(chromosomeN);
			dna.add(chrom);
			chromosomeN = new ArrayList<>();
		}
		
		ploidy = chromosomeN.size() + 1;
	}

	public void display(){
		Chromosome c;
		for(int i = 0; i < dna.size(); i++){
			c = dna.get(i);
			System.out.printf("Chromosome %d\n", i + 1);
			c.displayNBases();
			c.setCodons();
			/*System.out.println("\nCodons");
			c.displayCodons();*/
			System.out.println("\nAlleles\n");
			c.setAlleles();
			c.displayAlleles();
			c.setGenes();
		}
	}
	
	public void set(){
		for(Chromosome c : dna){
			c.setCodons();
			c.setAlleles();
			c.setGenes();
		}
	}
	
	public ArrayList<ArrayList<ArrayList<Character>>> getSequence(){
		ArrayList<ArrayList<ArrayList<Character>>> sequence = new ArrayList<>();
		for(Chromosome c : dna){
			sequence.add(c.getNBases());
		}
		return sequence;
	}
	
	ArrayList<Chromosome> dna = new ArrayList<>(); 
	
	public boolean isFemale(){ //completely arbitrary and inaccurate
		return false;
	}

	public void mutateAdd() {
		for(Chromosome c : dna){
			c.mutateAdd();
		}
	}
	
	public void mutateAdd(int prob){
		for(Chromosome c : dna){
			c.mutateAdd(prob);
		}
	}

	public void mutateDelete() {
		for(Chromosome c : dna){
			c.mutateDelete();
		}
	}
	
	public void mutateDelete(int prob) {
		for(Chromosome c : dna){
			c.mutateDelete(prob);
		}
	}
	
	public int getGeneNo(){
		int n = 0;
		for(Chromosome c : dna){
			n += c.getGeneNo();
		}
		return n;
	}
	
	public int size(){
		return dna.size();
	}
	
	public int ploidy(){
		return ploidy;
	}

	public void addChromosome() {
		Chromosome c = new Chromosome(ploidy);
		c.setRand(CHROMSIZE, GENEPERCHROM);		
		dna.add(c);
	}

	public void addPloid() {
		ploidy++;
		System.out.println("Adding ploid");
		for(Chromosome c : dna){
			c.setPloidy(ploidy);
		}
	}
	
	//private void getCharacteristics() //create phylogeneticTraits object
}

class Chromosome{
	public Chromosome(int ploid){
		if(ploid > 0){
			ploidy = ploid;
		} else {
			throw new IllegalArgumentException(
					"Ploidy needs to be 1 or higher");
		}
		/*insertStart(1);
		insertStart(2);
		insertRandomCodons(1, 2);
		insertRandomCodons(2, 2);
		insertStop(1);
		insertStop(2);*/
		//setGenes();
	}
	
	public Chromosome(ArrayList<ArrayList<Nucleotide>> chromosomes){
		ploidy = chromosomes.size() + 1;
		this.chromosomes = chromosomes;
	}
	
	private int ploidy;
	private ArrayList<ArrayList<Nucleotide>> chromosomes = new ArrayList<>(); //ploidy
	private ArrayList<ArrayList<Codon>> codons = new ArrayList<>();
	private ArrayList<ArrayList<Allele>> alleles = new ArrayList<>();
	private ArrayList<Gene> genes = new ArrayList<>();

	public int getPloidy() {
		return ploidy;
	}

	public void setPloidy(int ploidy) {
		this.ploidy = ploidy;
		/*ArrayList<Allele> chromatid = new ArrayList<>();
		Allele copy;
		System.out.println("before: " + alleles.size());
		for(Allele a : alleles.get(alleles.size() - 1)){
			copy = new Allele(a);
			chromatid.add(copy);
		}*/
		
		Nucleotide copy;
		ArrayList<Nucleotide> chromatid = new ArrayList<>();
		
		for(Nucleotide n : chromosomes.get(chromosomes.size() - 1)){
			copy = new Nucleotide(n.getNBase());
			chromatid.add(copy);
		}
		chromosomes.add(chromatid);
		
		//alleles.add(chromatid);
		setAlleles();
		setGenes();
		//System.out.println("after:  " + alleles.size());
	}

	public ArrayList<ArrayList<Nucleotide>> getNucleotides() {
		return chromosomes;
	}
	
	public ArrayList<ArrayList<Character>> getNBases(){
		ArrayList<ArrayList<Character>> nBases = new ArrayList <>();
		ArrayList<Character> chromosome = new ArrayList<>();
		
		for(ArrayList<Nucleotide> chromo : chromosomes){
			for(Nucleotide nBase : chromo){
				chromosome.add(nBase.getNBase());
			}
			nBases.add(chromosome);
			chromosome = new ArrayList<>();
		}
		return nBases;
	}
	
	public int getGeneNo(){
		return genes.size();
	}

	public void setNucleotides(ArrayList<ArrayList<Nucleotide>> chromosomes) {
		this.chromosomes = chromosomes;
	}

	private void insertStop(int chromo) { //TAA TAG TGA
		//maybe put error check for chromo		
		Random rand = new Random();
		int i;
		Nucleotide[] n = new Nucleotide[3];
		i = rand.nextInt(3);
		switch(i){
		case 0:
			n[0] = new Nucleotide('T');
			n[1] = new Nucleotide('A');
			n[2] = new Nucleotide('A');
			break;
		case 1:
			n[0] = new Nucleotide('T');
			n[1] = new Nucleotide('A');
			n[2] = new Nucleotide('G');
			break;
		case 2:
			n[0] = new Nucleotide('T');
			n[1] = new Nucleotide('G');
			n[2] = new Nucleotide('A');
			break;
		}
		
		Codon c = new Codon(n);
		codons.get(chromo).add(c);
	}

	private void insertStart(int chromo){
		ArrayList<Nucleotide> start = new ArrayList<>();
		start.addAll(Arrays.asList(stringToN("TATAAAAATG")));
		chromosomes.get(chromo).addAll(start);
	}
	
	private Nucleotide[] stringToN(String s){
		Nucleotide[] ns = new Nucleotide[s.length()];
		Nucleotide n;
		for(int i = 0; i < ns.length; i++){
			n = new Nucleotide(s.charAt(i));
			ns[i] = n;
		}
		return ns;
	}
	
	private void insertRandomNucleotides(int chromo, int quant) {
		ArrayList<Nucleotide> n;
		
		for(int i = 0; i < quant * 3; i++){ //3 bases per codon
			if(chromosomes.isEmpty()){
				n = new ArrayList<>();
				n.add(getRandomNucleotide());
				chromosomes.add(n);
			} else {
				chromosomes.get(chromo).add(getRandomNucleotide());
			}
		}
	}
	
	private Nucleotide getRandomNucleotide(){
		Nucleotide n;
		Random rand = new Random();
		int r;
		r = rand.nextInt(4);
		
		switch (r){
		case 0:
			n = new Nucleotide('A');
			return n;
		case 1:
			n = new Nucleotide('C');
			return n;
		case 2:
			n = new Nucleotide('T');
			return n;
		case 3:
			n = new Nucleotide('G');
			return n;
		default:
			return getRandomNucleotide();
		}
		

	}
	
	public void setRand(int chromSize, int genes) {
		Nucleotide n;
		Random rand = new Random();
		ArrayList<Nucleotide> nAL = new ArrayList<>();
		for(int i = 0; i < ploidy; i++){
			chromosomes.add(nAL);
			for(int j = 0; j < chromSize; j++){
				if(rand.nextInt(chromSize / genes) == 0){ //expected value of the number of genes as argument per chromosome
					insertStart(i);
				} else {
					n = new Nucleotide();
					chromosomes.get(i).add(n);
				}
			}
			nAL = new ArrayList<>();
		}
		setCodons();
		setAlleles();
		setGenes();
	}

	public void setCodons(){ //transcription
		codons = new ArrayList<>();
		Nucleotide[] codon = new Nucleotide[3];
		ArrayList<Nucleotide> chromo;
		String start = "TATAAAA";
		Codon c;
		boolean foundStart = false;
		int j = 0;

		ArrayList<Codon> cods = new ArrayList<>();
		
		//for(ArrayList<Nucleotide> chromo : chromosomes){
		for(int n = 0; n < chromosomes.size(); n++){
			chromo = chromosomes.get(n);
			for(int i = 0; i < chromo.size(); i++){ //TATAAAAA, cannot be programmed otherwise due to its arbitrary nature
				
				if(!foundStart){
					if(chromo.get(i).getNBase() == start.charAt(j)){
						j++;
						if(j == start.length()){
							foundStart = true;
							j = 0;
							i++; //skip last letter of start
						}
					} else {
						j = 0;
					}
				}
				//System.out.printf("\nchromo.get(i) = %c, start.charAt(j) = %c, i = %d j = %d", chromo.get(i).getNBase(), start.charAt(j), i, j);

								
				if(foundStart && i + 3 <= chromo.size()){
					//System.out.println("Found start");
					codon[0] = chromo.get(i);
					codon[1] = chromo.get(i + 1);
					codon[2] = chromo.get(i + 2);
					c = new Codon(codon);
					//System.out.printf("\nChromosome %d: %s, %c %c %c", n + 1, c.getName(), codon[0].getNBase(), codon[1].getNBase(), codon[2].getNBase());
					codon = new Nucleotide[3];
					if(c.getName() == "Stop"){
						foundStart = false;
					}
					cods.add(c);
					i += 2;
				}
			}
			codons.add(cods);
			cods = new ArrayList<>();
			foundStart = false;
		}
		
	}
	
	public void setAlleles(){
		alleles = new ArrayList<>();
		ArrayList<Codon> allele = new ArrayList<>();
		ArrayList<Allele> allelesInChrom = new ArrayList<>();
		Codon codon;
		Allele a;
		boolean foundStart = false;
		
		for(ArrayList<Codon> chromo : codons){ //chromosome in codon form
			for(int i = 0; i < chromo.size(); i ++){
				codon = chromo.get(i);
				
				if(codon.getName().equals("Methionine") && foundStart == false){
					foundStart = true;
					i++;
					if(i == chromo.size())
						break;
					codon = chromo.get(i);
				}
				
				
				if(foundStart){
					if(codon.getName().equals("Stop")){
						foundStart = false;
						if(allele.size() != 0){
							a = new Allele(allele);
							allelesInChrom.add(a);
							allele = new ArrayList<>();
						}
					} else {
						allele.add(codon);
						if(i == chromo.size() - 1){
							a = new Allele(allele);
							allelesInChrom.add(a);
							allele = new ArrayList<>();
						}
					}
				}
			}
			alleles.add(allelesInChrom);
			allelesInChrom = new ArrayList<>();
		}
	}
	
	public void setGenes(){
		ArrayList<Allele> gene = new ArrayList<>();
		Gene g;
		
		for(int j = 0; j < alleles.get(0).size(); j++){
			for(int i = 0; i < alleles.size(); i++){
				if(j < alleles.size()){
					gene = alleles.get(j);
					g = new Gene(gene);
					genes.add(g);
				}
			}
		}
	}
	
	public void displayNBases(){
		for(int i = 0; i < chromosomes.size(); i++){
			System.out.printf("\nChromatid %d\n", i + 1);
			for(Nucleotide n : chromosomes.get(i)){
				System.out.print(n.getNBase());
			}
		}
	}
	
	
	public void displayCodons(){
		for(int i = 0; i < codons.size(); i++){
			System.out.printf("\nChromatid %d of %d\n", i + 1, codons.size());
			for(Codon c: codons.get(i)){
				System.out.println(c.getName());
			}
		}
	}
	
	public void displayAlleles(){
		for(int i = 0; i < alleles.size(); i++){
			System.out.printf("Chromatid %d of %d\n", i + 1, alleles.size());
			for(Allele a: alleles.get(i)){
				System.out.printf("Size:%d\n", a.size());
				a.display();
				System.out.println();
			}
		}
	}
	
	public void mutateAdd(){
		Random rand = new Random();
		int r;
		for(ArrayList<Nucleotide> nal : chromosomes){
			for(int i = 0; i < nal.size(); i++){
				r = rand.nextInt(nal.size());
				if(r == 0){
					nal.add(i, getRandomNucleotide());
				}
			}
		}
	}
	
	public void mutateAdd(int prob){
		Random rand = new Random();
		int r;
		for(ArrayList<Nucleotide> nal : chromosomes){
			for(int i = 0; i < nal.size(); i++){
				r = rand.nextInt(prob);
				if(r == 0){
					nal.add(i, getRandomNucleotide());
				}
			}
		}
	}
	
	public void mutateDelete(){
		Random rand = new Random();
		int r;
		for(ArrayList<Nucleotide> nal : chromosomes){
			for(int i = 0; i < nal.size(); i++){
				r = rand.nextInt(nal.size());
				if(r == 0){
					nal.remove(i);
				}
			}
		}
	}
	
	public void mutateDelete(int prob){
		Random rand = new Random();
		int r;
		for(ArrayList<Nucleotide> nal : chromosomes){
			for(int i = 0; i < nal.size(); i++){
				r = rand.nextInt(prob);
				if(r == 0){
					nal.remove(i);
				}
			}
		}
	}
	
}

class Gene{
	
	public Gene(){
		throw new IllegalArgumentException(
				"Gene cannot be empty");
	}
	
	public Gene(ArrayList<Allele> a){
		for(Allele al : a){
			if(al == null){
				throw new IllegalArgumentException(
						"An allele is empty");
			}
		}
		
		alleles.addAll(a);
		ploidy = a.size() + 1;
	}
	
	private int ploidy;
	
	private ArrayList<Allele> alleles = new ArrayList<>();
	private ArrayList<String> definitions = new ArrayList<>();
}

class Allele{	
	public Allele(ArrayList<Codon> codon){ //only exons
		allele.addAll(codon);
	}
	
	public Allele(Allele a) {
		Nucleotide[] copyN = new Nucleotide[3];
		Codon copyC;
		for(Codon c : a.getCodons()){
			for(int i = 0; i < c.getCode().length; i++){
				copyN[i] = new Nucleotide(c.getCode()[i].getNBase());
			}
			copyC = new Codon(copyN);
			allele.add(copyC);
			copyN = new Nucleotide[3];
			copyC = new Codon(copyN);
		}
	}

	private ArrayList<Codon> allele = new ArrayList<>();
	private int dominance;
	
	public int size(){
		return allele.size();
	}
	
	public void addCodon(Codon c){
		allele.add(c);
	}
	
	public ArrayList<Codon> getCodons(){
		return allele;
	}
	
	public void display(){
		for(Codon c : allele){
			System.out.println(c.getName());
		}
	}
}

class Codon{	//or amino acid
	private static final int size = 3;
	
	public Codon(Nucleotide[] n){
		if(n.length == 0){
			throw new IllegalArgumentException(
					"Codon cannot be empty");
		}
		
		if(n.length != 3){
			throw new IllegalArgumentException(
					"Codon cannot be a size other than 3 neucleotides");
		}
		
		nucleotides = n;
	}
	
	public Codon(char c1, char c2, char c3){
		Nucleotide n1 = new Nucleotide(c1);
		Nucleotide n2 = new Nucleotide(c2);
		Nucleotide n3 = new Nucleotide(c3);
		Nucleotide[] n = {n1, n2, n3};
		nucleotides = n;
	}
	
	public Codon(Nucleotide n1, Nucleotide n2, Nucleotide n3){
		Nucleotide[] n = {n1, n2, n3};
		nucleotides = n;
	}
	
	private Nucleotide nucleotides[] = new Nucleotide[size]; 
	
	public Nucleotide[] getCode(){
		return nucleotides;
	}
	
	private static String names[][][] = { //not dealing with mRNA so Uracil is replaced with Thymine, 0:A 1:C 2:T 3:G
											{
												{"Lysine", "Asparagine", "Asparagine", "Lysine"},					  //AAX
												{"Threonine", "Threonine", "Threonine", "Threonine"},				  //ACX
												{"Isoleucine", "Isoleucine", "Isoleucine", "Methionine"},			  //ATX
												{"Argenine", "Serine", "Serine", "Arginine"} 						  //AGX
											},
											{
												{"Glutamine", "Histidine", "Histidine", "Glutamine"},				  //CAX
												{"Proline", "Proline", "Proline", "Proline"},						  //CCX
												{"Leucine", "Leucine", "Leucine", "Leucine"},						  //CTX
												{"Arginine", "Arginine", "Arginine", "Arginine"}					  //CGX
											},
											{
												{"Stop", "Tyrosine", "Tyrosine", "Stop"},							  //TAX
												{"Sernine", "Sernine", "Sernine", "Sernine"},						  //TCX
												{"Leucine", "Phenylalanine", "Phenylalanine", "Leucine"},			  //TTX
												{"Stop", "Cysteine", "Cysteine", "Tryptophan"}						  //TGX
											},
											{
												{"Glutamic acid", "Aspartic acid", "Aspartic acid", "Glutamic acid"}, //GAX
												{"Alanine", "Alanine", "Alanine", "Alanine"},						  //GCX
												{"Valine", "Valine", "Valine", "Valine"},							  //GTX
												{"Glycine", "Glycine", "Glycine", "Glycine"}						  //GGX
											}
										};
	public String getName(){
		int[] coord = new int [size];
		
		for(int i = 0; i < size; i++){
			switch(nucleotides[i].getNBase()){
			case 'A':
				coord[i] = 0;
				break;
			case 'C':
				coord[i] = 1;
				break;
			case 'T':
				coord[i] = 2;
				break;
			case 'G':
				coord[i] = 3;
				break;
			}
		}
		return names[coord[0]] [coord[1]] [coord[2]];
	}
}

class Nucleotide{
	
	public Nucleotide(char n){
		if(validNucleotide(n)){
			setNucleotide(n);
		} else {
			throw new IllegalArgumentException(
					"Invalid Codon");
		}
		if(nBase == '\u0000'){
			throw new IllegalArgumentException(
					"Neucleotide cannot be null");
		}
	}
	
	public Nucleotide(){
		setRand();
	}
	
	private boolean validNucleotide(char n) {
		if (n == 'A' || n == 'C' || n == 'T' || n == 'G'){
			return true;
		} else {
			return false;
		}
	}
	
	public void setNucleotide(char n) {
		nBase = n;
	}
	
	public char getNBase(){
		return nBase;
	}
	
	private Random rand = new Random();
	
	public void setRand(){
		int r;
		r = rand.nextInt(4);
		switch(r){
		case 0:
			nBase = 'A';
			break;
		case 1:
			nBase = 'C';
			break;
		case 2:
			nBase = 'T';
			break;
		case 3:
			nBase = 'G';
			break;
		}
	}
	
	public boolean equals(char c){
		if(c == nBase){
			return true;
		} else {
			return false;
		}
		
	}

	private char nBase;
	
}