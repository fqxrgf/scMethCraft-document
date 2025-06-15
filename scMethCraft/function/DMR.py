def plot_DMR(adata,focus_type):
    result = adata.uns['rank_genes_groups']

    data = pd.DataFrame({
        'log2_fold_change': result["logfoldchanges"][focus_type],
        'p_value': result["pvals"][focus_type]
    })


    data['neg_log10_pvalue'] = -np.log10(data['p_value'])


    p_value_threshold = 0.01
    log2fc_threshold = 1


    data['significance'] = 'Not significant'
    data.loc[:, 'significance'] = 'None'
    data.loc[(data['p_value'] < p_value_threshold) & 
             (data['log2_fold_change'] <= -log2fc_threshold), 'significance'] = 'Hypo DMR'
    data.loc[(data['p_value'] < p_value_threshold) & 
             (data['log2_fold_change'] >= log2fc_threshold), 'significance'] = 'Hyper DMR'


    # 绘制火山图
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=data,
        x='log2_fold_change',
        y='neg_log10_pvalue',
        hue='significance',
        palette={
            'None': '#a3c5df',
            'Low': '#abc996',
            'High': '#eb95c6'
        },
        alpha=0.8,           
        s=60,                
        edgecolor='white',   
        linewidth=0.3        
    )



    plt.axhline(-np.log10(p_value_threshold), color='black', linestyle='--', linewidth=1)
    plt.axvline(-log2fc_threshold, color='black', linestyle='--', linewidth=1)
    plt.axvline(log2fc_threshold, color='black', linestyle='--', linewidth=1)

    plt.title('Volcano Plot', fontsize=16)
    plt.xlabel('log2(Fold Change)', fontsize=14)
    plt.ylabel('-log10(p-value)', fontsize=14)
    plt.legend(title='Significance', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()
