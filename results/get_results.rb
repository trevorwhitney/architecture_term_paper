Dir[File.dirname(__FILE__) + '/*.rb'].each { |f| require f }

results = Results.new("../O0/final_results.csv", "../O1/final_results.csv", 
  "../O2/final_results.csv", "../O3/final_results.csv", 
  "../O3_funroll/final_results.csv")

results.generate_results