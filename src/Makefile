build:
	javac *.java

run: build
	java Search params/degrees_10.params > ./summaries/degrees.out
	mv ./pointscatter_summary.txt ./summaries/degrees.txt
	java Search params/thetamag_10.params > ./summaries/thetamag.out
	mv ./pointscatter_summary.txt ./summaries/thetamag.txt
	java Search params/xycart_10.params > ./summaries/xycart.out
	mv ./pointscatter_summary.txt ./summaries/xycart.txt

degrees: build
	java Search params/degrees_10.params

thetamag: build
	java Search params/thetamag_10.params

xycart: build
	java Search params/xycart_10.params

experiment: build
	java Search params/degrees_10.params > ./summaries/degrees_10.out
	mv ./pointscatter_summary.txt ./summaries/degrees_10.txt
	mv ./best_fit.csv ./summaries/degrees_10.csv
	mv ./average_stats.csv ./summaries/degrees_10_averages.csv
	java Search params/thetamag_10.params > ./summaries/thetamag_10.out
	mv ./pointscatter_summary.txt ./summaries/thetamag_10.txt
	mv ./best_fit.csv ./summaries/thetamag_10.csv
	mv ./average_stats.csv ./summaries/thetamag_10_averages.csv
	java Search params/xycart_10.params > ./summaries/xycart_10.out
	mv ./pointscatter_summary.txt ./summaries/xycart_10.txt
	mv ./best_fit.csv ./summaries/xycart_10.csv
	mv ./average_stats.csv ./summaries/xycart_10_averages.csv
	java Search params/degrees_25.params > ./summaries/degrees_25.out
	mv ./pointscatter_summary.txt ./summaries/degrees_25.txt
	mv ./best_fit.csv ./summaries/degrees_25.csv
	mv ./average_stats.csv ./summaries/degrees_25_averages.csv
	java Search params/thetamag_25.params > ./summaries/thetamag_25.out
	mv ./pointscatter_summary.txt ./summaries/thetamag_25.txt
	mv ./best_fit.csv ./summaries/thetamag_25.csv
	mv ./average_stats.csv ./summaries/thetamag_25_averages.csv
	java Search params/xycart_25.params > ./summaries/xycart_25.out
	mv ./pointscatter_summary.txt ./summaries/xycart_25.txt
	mv ./best_fit.csv ./summaries/xycart_25.csv
	mv ./average_stats.csv ./summaries/xycart_25_averages.csv
	java Search params/degrees_100.params > ./summaries/degrees_100.out
	mv ./pointscatter_summary.txt ./summaries/degrees_100.txt
	mv ./best_fit.csv ./summaries/degrees_100.csv
	mv ./average_stats.csv ./summaries/degrees_100_averages.csv
	java Search params/thetamag_100.params > ./summaries/thetamag_100.out
	mv ./pointscatter_summary.txt ./summaries/thetamag_100.txt
	mv ./best_fit.csv ./summaries/thetamag_100.csv
	mv ./average_stats.csv ./summaries/thetamag_100_averages.csv
	java Search params/xycart_100.params > ./summaries/xycart_100.out
	mv ./pointscatter_summary.txt ./summaries/xycart_100.txt
	mv ./best_fit.csv ./summaries/xycart_100.csv
	mv ./average_stats.csv ./summaries/xycart_100_averages.csv

clean:
	rm -f *.class
