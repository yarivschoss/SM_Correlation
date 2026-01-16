function main()
    % MAIN – entry point for the transformer–household mapping project

    % Add all subfolders to the MATLAB path
    setup_paths();

    % Load configuration for data collection / simulation
    cfg = get_data_config();

    % Step 1: Load or generate raw DATA
    rawData = load_or_generate_data(cfg);

    % Step 2: Pre-processing of raw data
    data = preprocess_data(rawData, cfg);

    plot_raw_vs_clean(rawData, data, "Transformer1", [datetime(2025,8,13) datetime(2025,8,14)]);

    % Step 3: Run energy-balance-based algorithm
    results.energyBalance = run_energy_balance(data, cfg);

    % Step 4: Run optimization (Gradient / Stochastic Gradient Descent)
    results.optimization = run_optimization(data, cfg);

    % Step 5: Run E-CPC algorithm for event-based connectivity
    results.ecpc = run_ecpc(data, cfg);

    % Step 6: Evaluate performance and display results
    evaluate_results(results, cfg);
end
