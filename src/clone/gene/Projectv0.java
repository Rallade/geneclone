package clone.gene;
import java.util.ArrayList;
import java.util.Random;

public class Projectv0 {
	private Random rand = new Random();
	private ArrayList<Entity> generation = new ArrayList<>();
	private ArrayList<ArrayList<Entity>> entities = new ArrayList<>();
	
	private final int MUTFREQ = 20; // mutation frequency, the lower, the more frequent

	public void start() {
		Entity LUCA = new Entity(1, 1); //last universal common ancestor
		generation.add(LUCA);
		entities.add(generation);
	}

	public void procreate() {
		Entity d1; // daughter
		Entity d2;
		generation = new ArrayList<>();
		for (Entity e : entities.get(entities.size() - 1)) { // last generation
			if (!e.isSexuallyReproducing()) {
				d1 = new Entity(e.getSequence());
				d2 = new Entity(e.getSequence());
				d1.mutate((generation.size() + 1) * MUTFREQ); // with each generation, entities become better at not mutating
				d2.mutate((generation.size() + 1) * MUTFREQ);
				if(d1.isViable()){
					generation.add(d1);
				}
				if(d2.isViable()){
					generation.add(d2);
				}
			} else {
				// mate(e, findMate(e));
			}
		}
		entities.add(generation);
	}

	private Entity findMate(Entity e) {
		return null;
	}

	public Entity mate(Entity parent1, Entity parent2) {
		return null;
	}
	
	public void display(){
		for(int i = 0; i < entities.size(); i++){
			System.out.println("Generation " + (i + 1) + ":");
			for(Entity e : entities.get(i)){
				e.display();
			}
		}
	}

	public static void main(String[] args) {
		Projectv0 p = new Projectv0();
		p.start();
		for(int i = 0; i < 4; i++){
			p.procreate();
		}
		p.display();
		
		//DNA dna = new DNA(1, 1);
		//dna.addChromosome();
		//dna.addPloid();
		//dna.display();
	}
}

class Entity {
	private Random rand = new Random();
	private final int CHROMSMOD = 2; // chromosome number increase frequency modifier (1/x)
	private final int PLOIDMOD = 3; //ploidy number increase frequency modifier (1/x)

	public Entity(int chromosomes, int homologuous) {
		if (chromosomes == 0 || homologuous == 0) {
			throw new IllegalArgumentException("Cannot have 0 chromosomes");
		}
		dna = new DNA(chromosomes, homologuous);
		this.chromosomes = chromosomes;
		this.homologuous = homologuous;
		setReproduction();
	}

	public boolean isViable() {
		if(dna.getGeneNo() > 0){
			return true;
		} else {
			return false;
		}
	}

	public void mutate() {
		if(rand.nextBoolean()){
			dna.mutateAdd();
		}
		if(rand.nextBoolean()){
			dna.mutateDelete();
		}
		set();
	}
	
	public void mutate(int freq) {
		if(rand.nextBoolean()){
			dna.mutateAdd(freq);
		}
		if(rand.nextBoolean()){
			dna.mutateDelete(freq);
		}
		if(rand.nextInt(CHROMSMOD * freq * dna.size()) == 0){ //less and less frequent
			dna.addChromosome();
		}
		if(rand.nextInt(PLOIDMOD * freq * dna.ploidy()) == 0){
			dna.addPloid();
		}
		set();
	}

	public void display() {
		dna.display();
		
	}

	public Entity(ArrayList<ArrayList<ArrayList<Character>>> sequence) {
		dna = new DNA(sequence);
	}

	private void setReproduction() {
		if (homologuous == 1 && chromosomes == 1) { // prokaryotes divide
			sexualReproduction = false;
		} else {
			sexualReproduction = true;
		}
	}

	public boolean isSexuallyReproducing() {
		return sexualReproduction;
	}

	public void setSexualReproduction(boolean sexualReproduction) {
		this.sexualReproduction = sexualReproduction;
	}

	public int getChromosomes() {
		return chromosomes;
	}

	public int getHomologuous() {
		return homologuous;
	}
	
	public void set(){
		dna.set();
	}

	public ArrayList<ArrayList<ArrayList<Character>>> getSequence() {
		return dna.getSequence();
	}

	private DNA dna;
	private boolean sexualReproduction;
	private int chromosomes;
	private int homologuous;

}