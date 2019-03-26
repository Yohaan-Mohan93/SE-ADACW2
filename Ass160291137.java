//Name: Yohaan Mohan
//SRN: 160291137
import java.util.HashSet;
import java.util.ArrayList;
import java.io.*;
import java.util.Scanner;
import java.util.*;

public class Ass160291137 {
	
	static int N= 1000;
	static double [][] edges = new double[N][N];
  static double [][] weightedEdges = new double[N][N];
  static double [][] distanceEdges = new double[N][N];
	static TreeMap <Integer,String> cityNames = new TreeMap<Integer,String>();
  static TreeMap <Integer,String> distances = new TreeMap<Integer,String>();
	
	public static class graph
	{
		double [] [] adj;
    //Provided by UOL
		graph (double [] [] a)
		{
			adj= new double [a.length][a.length];
			for (int i=0;i<a.length;i++)
				for (int j=0;j<a.length;j++) adj[i][j]=a[i][j];
		}
		//Provided by UOL
    @SuppressWarnings("unchecked")
    int findSmallest(HashMap <Integer,Double> t)
    {
      Object [] things= t.keySet().toArray();
      double val=t.get(things[0]);
      int least=(int) things[0];
      Set <Integer> k = t.keySet();
      
      for (Integer i : k)
      {
        if (t.get(i) < val)
        {
          least=i;
          val=t.get(i);
        }
      }
      return least;
    }
    //Provided by UOL, finished by me
    @SuppressWarnings("unchecked")
    public ArrayList <Integer> dijkstra (int start, int end)
    {
      int N= 1000;
      HashMap <Integer,Double> Q = new HashMap <Integer,Double>();
      ArrayList <Integer> [] paths = new ArrayList [N];
    
      for (int i=0;i<N;i++)
      {
        Q.put(i, Double.POSITIVE_INFINITY);
        paths[i]=new ArrayList <Integer>();
        paths[i].add(start);
      }
      HashSet <Integer> S= new HashSet();
      S.add(start);
      Q.put(start,0.0);
      while (!Q.isEmpty())
      {   
        int v = findSmallest(Q);
        if (v==end && (Q.get(v) != Double.POSITIVE_INFINITY))
          return paths[end];
        double w = Q.get(v);  
        S.add(v);
        for(int u: neighbours(v))
        if (!(S.contains(u)))
        {
          double w1= adj[v][u] + w;
          if (w1 < Q.get(u))
          {
            Q.put(u,w1);
            paths[u]= addToEnd(u, paths[v]);
          }
        }
      Q.remove(v);
      }
      return new ArrayList <Integer> ();
    }
    //Provided by UOL
		public HashSet <Integer> neighbours(int v)
		{
			HashSet <Integer> h = new HashSet <Integer> ();
			for (int i=0;i<adj.length;i++) if (adj[v][i]!=0) h.add(i);
			return h;
		}
		//Provided by UOL
		public HashSet <Integer> vertices()
		{
			HashSet <Integer> h = new HashSet <Integer>();
			for (int i=0;i<adj.length;i++) h.add(i);
			return h;
		}
		
    //Provided by UOL
		@SuppressWarnings("unchecked")
		ArrayList <Integer> addToEnd (int i, ArrayList <Integer> path)
		// returns a new path with i at the end of path
		{
			ArrayList <Integer> k; k=(ArrayList<Integer>)path.clone(); k.add(i);
			return k;
		}
		
    //Provided by UOL
		@SuppressWarnings("unchecked")
		public HashSet <ArrayList <Integer>> shortestPaths1(HashSet<ArrayList <Integer>> sofar, 
                                                         HashSet <Integer> visited, int end)
		{
			HashSet <ArrayList <Integer>> more = new HashSet <ArrayList<Integer>>();
			HashSet <ArrayList <Integer>> result = new HashSet <ArrayList<Integer>>();
			HashSet <Integer> newVisited = (HashSet <Integer>)visited.clone();
			boolean done = false;
			boolean carryon = false;
			for (ArrayList <Integer> p: sofar)
			{
				for (Integer z: neighbours(p.get(p.size()-1)))
				{
					if (!visited.contains(z))
					{
						carryon=true; newVisited.add(z);
						if (z==end) {
							done=true;
							result.add(addToEnd(z,p));
						}
						else
							more.add(addToEnd(z,p));
					}
				}
			}
			if (done) return result; else
				if (carryon)
					return shortestPaths1(more,newVisited,end);
				else
					return new HashSet <ArrayList <Integer>>();
		}
		
    //Provided by UOL
		public HashSet <ArrayList <Integer>> shortestPaths( int first, int end)
		{
			HashSet <ArrayList <Integer>> sofar = new HashSet <ArrayList<Integer>>();
			HashSet <Integer> visited = new	HashSet<Integer>();
			ArrayList <Integer> starting = new ArrayList<Integer>();
			starting.add(first);
			sofar.add(starting);
			if (first==end)
			return sofar; 
			visited.add(first);
			return shortestPaths1(sofar,visited,end);
		}
  }
  
  //Provided by UOL
  static double realDistance(double lat1, double lon1, double lat2, double lon2)
  {
    int R = 6371;
    // km (change this constant to get miles)
    double dLat = (lat2-lat1) * Math.PI / 180;
    double dLon = (lon2-lon1) * Math.PI / 180;
    double a = Math.sin(dLat/2) * Math.sin(dLat/2) +
    Math.cos(lat1 * Math.PI / 180 ) * Math.cos(lat2 * Math.PI / 180 )* Math.sin(dLon/2) * Math.sin(dLon/2);
    double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a)); double d = R * c;
    return d;
  }
  
  //Provided by UOL  
  @SuppressWarnings("unchecked")
  static ArrayList<Integer> firstElement (HashSet <ArrayList <Integer>> s)
  {
    return ( ArrayList<Integer>)s.toArray()[0];
  }
  //Provided by UOL
	static ArrayList<String> convert(ArrayList<Integer> m)
	{
		ArrayList<String> z= new ArrayList<String>();
		for (Integer i:m)
			z.add(cityNames.get(i));
		return z;
	}
	//Provided by UOL
	static HashSet<ArrayList<String>> convert (HashSet<ArrayList<Integer>> paths)
	{
		HashSet <ArrayList <String>> k= new HashSet<ArrayList<String>>();
		for (ArrayList <Integer> p:paths)
			k.add(convert(p));
		return k;
	}
	
	public static void main(String[] args) throws Exception {
    int q1Ans, question1City1, question1City2, highest, question2city1, question2city2;
    int  question5City, q5High, question6City1, question6City2, question7City1, question7City2;
    HashSet<ArrayList<Integer>> question4;
    Iterator<ArrayList<Integer>> iterator;
    Vector<Integer> pathSizes, question5Ans;
    TreeMap<Integer, Integer> question5;
    ArrayList<Integer> q6Ans, q7Ans;
    
    question1City1 = question1City2 = highest = question2city1 = question2city2 = question5City = 0;
    q5High = question6City1 = question6City2 = question7City1 = question7City2 = 0;
		System.out.println("Name: Yohaan Mohan");
		System.out.println("Student ID: 160291137");
		long startTime = System.nanoTime();
		
    for(int i=0;i<N;i++)
			for(int j=0;j<N;j++)
				edges[i][j]=0.0;
      
		Scanner s = new Scanner(new FileReader(args[2]));
		String z =s.nextLine();
		while (s.hasNext())
		{
			z =s.nextLine();
			String[] results = z.split(",");
			edges[Integer.parseInt(results[0])][Integer.parseInt(results[1])]=1.0;
			edges[Integer.parseInt(results[1])]	[Integer.parseInt(results[0])]=1.0;
		}
		
		s = new Scanner(new FileReader(args[0]));
		z =s.nextLine();
		while (s.hasNext())
		{
			z =s.nextLine();
			String[] results = z.split(",");
      if(results[1].equals("Athens")) {
				question1City1 = Integer.parseInt(results[0]);
			}
			if(results[1].equals("Tehran")) {
				question1City2 = Integer.parseInt(results[0]);
			}
      if(results[1].equals("Toronto")){
        question5City = Integer.parseInt(results[0]);
      }
      if(results[1].equals("Rome")){
        question6City1 = Integer.parseInt(results[0]);
      }
      if(results[1].equals("New York")){
        question6City2 = Integer.parseInt(results[0]);
      }
      if(results[1].equals("Lisbon")){
        question7City1 = Integer.parseInt(results[0]);
      }
      if(results[1].equals("Manila")){
        question7City2 = Integer.parseInt(results[0]);
      }
			cityNames.put(Integer.parseInt(results[0]),results[1]);
		}
		s.close();
		graph G = new graph(edges);
        
		//Work needed for question 1
    q1Ans = G.shortestPaths(question1City1,question1City2).size();
    
    //Work needed for question 2 and 3
    Vector<Integer> tempV = new Vector<Integer>();
    for(Map.Entry<Integer,String> i: cityNames.entrySet()){
      Integer key = i.getKey();
      String value = i.getValue();
      tempV.add(key);
      for(Map.Entry<Integer,String> j: cityNames.entrySet()){
        if(tempV.contains(j.getKey()))
          continue;
        int temp = G.shortestPaths(i.getKey(),j.getKey()).size();
        if(highest < temp) {
          highest = temp;
          question2city1 = i.getKey();
          question2city2 = j.getKey();
        }
      }
    }
    
    //Work neeeded for question 4
    question4 = G.shortestPaths(question2city1,question2city2);
    pathSizes = new Vector<Integer>();
    iterator = question4.iterator();
    while(iterator.hasNext()){
        pathSizes.add(firstElement(question4).size());
        iterator.next();
        iterator.remove();
    }
    
    //Work needed for question 5
    question5 = new TreeMap<Integer, Integer>();
    question5Ans = new Vector<Integer>();
    for(Map.Entry<Integer,String> i: cityNames.entrySet()){
      Integer key = i.getKey();
      String value = i.getValue();
      question5.put(key,G.shortestPaths(question5City,key).size());
    }
    
    for(Map.Entry<Integer,Integer> i: question5.entrySet()){
      if(q5High < i.getValue())
      {
        q5High = i.getValue(); 
      }
    }
    
    for(Map.Entry<Integer,Integer> i: question5.entrySet()){
      if(q5High == i.getValue())
      {
        question5Ans.add(i.getKey()); 
      }
    }
    
    //Create a weighted graph
    for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      weightedEdges[i][j]=0.0;
    
    s = new Scanner(new FileReader(args[2]));
    while (s.hasNext()){
      z =s.nextLine();
      String[] results = z.split(",");
      weightedEdges[Integer.parseInt(results[0])]
      [Integer.parseInt(results[1])]=Double.parseDouble(results[2]);
    }
    graph dijGraph = new graph(weightedEdges);
    
    //Work needed for question 6
    double q6sum = 0;
    q6Ans = new ArrayList<Integer>();
    q6Ans = dijGraph.dijkstra(question6City1, question6City2);
    if(q6Ans.size() == 2){
      HashSet<Integer> city1Neighbors = new HashSet<Integer>();
      HashSet<Integer> city2Neighbors = new HashSet<Integer>();
      TreeMap<Double, Integer> weightSums = new TreeMap<Double, Integer>();
      ArrayList<Integer> common = new ArrayList<Integer>();
      double lowest = 0;
      int c1 = question6City1;
      int c2 = question6City2;
      
      city1Neighbors = G.neighbours(q6Ans.get(0));
      city2Neighbors = G.neighbours(q6Ans.get(1));
      for(int i : city1Neighbors){
        if(city2Neighbors.contains(i)){
          common.add(i);
        }
      }
      for(int i = 1; i < common.size(); i++){
       if(weightedEdges[common.get(i)][c2] == 0 || weightedEdges[c1][common.get(i)] == 0)
         continue;
       else 
        weightSums.put((weightedEdges[c1][common.get(i)] + weightedEdges[common.get(i)][c2]),
                        common.get(i));
      }
      q6sum = weightSums.firstEntry().getKey();
    }
    
    //Create a distance graph
    for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      distanceEdges[i][j]=0.0;
    
    s = new Scanner(new FileReader(args[1]));
    while (s.hasNext()){
      z =s.nextLine();
      String[] results = z.split(",");
      distances.put(Integer.parseInt(results[0]),new String(results[1] +","+results[2]));
    }

    for(int i=0;i<N;i++){
      if(distances.containsKey(i)){
        for(int j=0;j<N;j++){
          if(distances.containsKey(j)){
            String c1 = distances.get(i);
            String c2 = distances.get(j);
            String[] c1LonLat = c1.split(",");
            String[] c2LonLat = c2.split(",");
            double c1Lon = Double.parseDouble(c1LonLat[0]);
            double c1Lat = Double.parseDouble(c1LonLat[1]);
            double c2Lon = Double.parseDouble(c2LonLat[0]);
            double c2Lat = Double.parseDouble(c2LonLat[1]);
            
            distanceEdges[i][j] = realDistance(c1Lon,c1Lat,c2Lon,c2Lat);
          }
          else
            continue;
        } 
      }
      else
        continue;
    }
    graph distGraph = new graph(distanceEdges);
    
    //Work needed for question 7
    double q7sum = 0;
    q7Ans = new ArrayList<Integer>();
    q7Ans = distGraph.dijkstra(question7City1,question7City2);
    if(q7Ans.size() == 2){
      HashSet<Integer> city1Neighbors = new HashSet<Integer>();
      TreeMap<Double, Integer> weightSums = new TreeMap<Double, Integer>();
      double lowest = 0;
      int c1 = question7City1;
      int c2 = question7City2;
      
      city1Neighbors = G.neighbours(q7Ans.get(0));
      for(int i : city1Neighbors){
        ArrayList<Integer> temp = new ArrayList<Integer>();
        temp = distGraph.dijkstra(i,question7City2);
        if(temp.size() == 0)
          continue;
        weightSums.put(distanceEdges[c1][temp.get(0)] + distanceEdges[temp.get(0)][temp.get(1)],i);
      }
      q7sum = weightSums.firstEntry().getKey();
    }
    
    //Question 1
    System.out.println("Question 1: " + q1Ans);
    
    //Question 2
    System.out.println("Question 2: " + question2city1 + " and " + question2city2);
    
    //Question 3
    System.out.println("Question 3: " + highest + " paths");
    
    //Question 4
    System.out.println("Question 4: " + pathSizes);
    
    //Question 5
    System.out.println("Question 5: " + question5Ans);
    
    //Question 6
    System.out.println("Question 6: " + q6sum);
    
    //Question 7
    System.out.println("Question 7: " + q7sum + " km");
        
		long endTime = System.nanoTime();
		long timeTaken = (endTime - startTime)/1000000;
		System.out.println("Execution Time: " + timeTaken + " milliseconds");
	}
}
