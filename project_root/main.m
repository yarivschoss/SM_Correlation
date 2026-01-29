function main()
    % MAIN – entry point for the transformer–household mapping project

    % Add all subfolders to the MATLAB path
    setup_paths();

    % Load configuration for data collection / simulation
    cfg = get_data_config();

    % Step 1: Run energy-balance-based algorithm
    results.energyBalance = run_energy_balance_ukf();

    % Step 2: Run optimization (Gradient / Stochastic Gradient Descent)
    results.optimization = run_optimization(data, cfg);

    % Step 3: Run E-CPC algorithm for event-based connectivity
    results.ecpc = run_ecpc(data, cfg);

    % Step 4: Evaluate performance and display results
    evaluate_results(results, cfg);
end
