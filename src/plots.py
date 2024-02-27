import plotly.graph_objects as go
import plotly
import numpy as np
import os
import itertools

def plot_coverage_data(html_graphs, arguments, chrom, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets, ref_end_values_phasesets, sufix):
    fig = go.Figure()
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=None, yaxis=None,
                               opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=None, yaxis=None,
                               opacity=0.7, color='steelblue')

    if arguments['unphased_reads_coverage_enable']:
        add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased', text=None,
                                   yaxis=None, opacity=0.7, color='olive')
    plots_add_markers_lines(fig)

    if arguments['phaseblocks_enable']:
        gaps_values = np.full(len(haplotype_1_values_phasesets), 'None')
        haplotype_1_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values)))
        haplotype_2_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values)))
        phaseblocks_positions = list(
            itertools.chain.from_iterable(zip(ref_start_values_phasesets, ref_end_values_phasesets, gaps_values)))

        add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values,
                                      haplotype_2_phaseblocks_values)

    plots_layout_settings(fig, chrom, arguments, ref_end_values[-1:][0], arguments['cut_threshold'])

    if arguments['pdf_enable']:
        print_chromosome_pdf(fig, chrom, arguments['out_dir_plots'])

    print_chromosome_html(fig, chrom + '_' + sufix, html_graphs, arguments['out_dir_plots'])
    html_graphs.write(
        "  <object data=\"" + chrom + '_' + sufix + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")


def plot_coverage_data_after_correction(html_graphs, arguments, chrom, ref_start_values, ref_end_values, haplotype_1_values, haplotype_2_values, unphased_reads_values, haplotype_1_values_phasesets, haplotype_2_values_phasesets, ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2, sufix):
    fig = go.Figure()
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_1_values, name='HP-1', text=None, yaxis=None,
                               opacity=0.7, color='firebrick')
    add_scatter_trace_coverage(fig, ref_start_values, haplotype_2_values, name='HP-2', text=None, yaxis=None,
                               opacity=0.7, color='steelblue')

    if arguments['unphased_reads_coverage_enable']:
        add_scatter_trace_coverage(fig, ref_start_values, unphased_reads_values, name='Unphased', text=None,
                                   yaxis=None, opacity=0.7, color='olive')
    plots_add_markers_lines(fig)

    if arguments['phaseblocks_enable']:
        gaps_values_hp1 = np.full(len(haplotype_1_values_phasesets), 'None')
        gaps_values_hp2 = np.full(len(haplotype_2_values_phasesets), 'None')
        haplotype_1_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_1_values_phasesets, haplotype_1_values_phasesets, gaps_values_hp1)))
        haplotype_2_phaseblocks_values = list(
            itertools.chain.from_iterable(zip(haplotype_2_values_phasesets, haplotype_2_values_phasesets, gaps_values_hp2)))
        phaseblocks_positions_hp1 = list(
            itertools.chain.from_iterable(zip(ref_start_values_phasesets_hp1, ref_end_values_phasesets_hp1, gaps_values_hp1)))
        phaseblocks_positions_hp2 = list(
            itertools.chain.from_iterable(zip(ref_start_values_phasesets_hp2, ref_end_values_phasesets_hp2, gaps_values_hp2)))

        add_scatter_trace_phaseblocks_seperate(fig, phaseblocks_positions_hp1, phaseblocks_positions_hp2, haplotype_1_phaseblocks_values,
                                      haplotype_2_phaseblocks_values)

    plots_layout_settings(fig, chrom, arguments, ref_end_values[-1:][0], arguments['cut_threshold'])

    if arguments['pdf_enable']:
        print_chromosome_pdf(fig, chrom, arguments['out_dir_plots'])

    print_chromosome_html(fig, chrom + '_' + sufix, html_graphs, arguments['out_dir_plots'])
    html_graphs.write(
        "  <object data=\"" + chrom + '_' + sufix + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")

def change_point_detection(data, start, ends, arguments, chrom, html_graphs, hp, color):
    import ruptures as rpt
    fig = go.Figure()
    #starts = [i for i in range(0, len(data), 50000)]
    add_scatter_trace_coverage(fig, start, data, name='HP-'+str(hp), text=None, yaxis=None,
                               opacity=0.7, color=color)

    data = np.array(data, dtype='int') #numpy.clip(data, a_min=0, a_max=1000)
    algo = rpt.Pelt(model="rbf", jump=25).fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i < len(data)]
    for i, point in enumerate(change_points):
        fig.add_vline(x=point*50000, y0=-10, y1=500, line_width=1, line_dash="dash",
                  line_color=color)

    plots_add_markers_lines(fig)
    plots_layout_settings(fig, chrom, arguments, ends[-1:][0], arguments['cut_threshold'])

    print_chromosome_html(fig, chrom + '_hp_'  + str(hp), html_graphs, arguments['out_dir_plots'])
    html_graphs.write("  <object data=\"" + chrom + '_hp_'  + str(hp)  + '.html' + "\" width=\"700\" height=\"420\"></object>" + "\n")


def print_chromosome_pdf(fig, chrom, coverage_plots_path):
    filename = f"{os.path.join(coverage_plots_path, chrom + '.pdf')}"
    plotly.io.write_image(fig,  filename, format='pdf')

def print_chromosome_html(fig, chrom, html_graphs, coverage_plots_path):
    fname = f"{os.path.join(coverage_plots_path, chrom + '.html')}"
    plotly.offline.plot(fig, filename=fname,auto_open=False)

def plots_add_markers_lines(fig):
    # style all the traces
    fig.update_traces(
        # hoverinfo="name+x+text+y",
        # line={"width": 0.5},
        marker={"size": 2},
        mode="markers",
        showlegend=True
    )

def add_scatter_trace_coverage(fig, x, y, name, text, yaxis, opacity, color, visibility=True):

    fig.add_trace(go.Scatter(
        #legendgroup="group1",  # this can be any string, not just "group"
        #legendgrouptitle_text="Coverage",
        x=x,
        y=y,
        name=name,
        text=text,
        #yaxis="y5",
        opacity=opacity,
        marker_color=color,
        visible=visibility,
    ))

def plots_layout_settings(fig, chrom, arguments, limit_x, limit_y):
    # Update axes
    fig.update_layout(
        xaxis=dict(
            #autorange=True,
            type="linear",
            showline=True,
            zeroline=True,
            linecolor = "dimgray",
            range=[0,limit_x*1.01]
        ),
        yaxis=dict(
            linecolor="dimgray",
            range=[0, limit_y+1],
            side="left",
            tickfont={"color": "dimgray"},
            tickmode="auto",
            ticks="outside",
            title="<b>Coverage</b> (mean depth)",
            titlefont={"color": "dimgray"},
            type="linear",
            showline=True,
            zeroline=True,
        ),
        yaxis2=dict(
            linecolor="dimgray",
            range=[1, arguments['cut_threshold'] + 5],
            side="right",
            tickfont={"color": "dimgray"},
            tickmode="auto",
            ticks="outside",
            title="<b>Copynumbers</b> (integers)",
            titlefont={"color": "dimgray"},
            type="linear",
            # showline=True,
            # zeroline=True,
            anchor="x",
            overlaying="y",
        )
    )
    # Update layout
    fig.update_layout(
        template="plotly_white",
        font_family="Times New Roman"
    )

    fig.update_layout(
        title=chrom,
    )
    #Legend
    fig.update_layout(legend=dict(
        orientation = 'h', xanchor = "center", x = 0.45, y= 1.2, #orientation = 'v', xanchor = "center", x = 1.08, y= .5
    ))
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    fig.update_xaxes(tick0=0.0, rangemode="nonnegative")

    fig.update_layout(legend={'itemsizing': 'constant'})

    fig.update_layout(font_family= "Times New Roman")

    fig.update_layout(
        title={
            'text': chrom + ' - ' +arguments['genome_name'],
            'y':0.96,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family = "Courier New",
        font_color = "dimgray",
        title_font_family = "Times New Roman",
        title_font_color = "red",
        legend_title_font_color = "green",
    )
    #Size
    fig.update_layout(
        width=680,
        height=400,
       )

def add_scatter_trace_phaseblocks_seperate(fig, phaseblocks_positions_hp1, phaseblocks_positions_hp2, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values):
    fig.add_trace(go.Scatter(
        #legendgroup="group3",  # this can be any string, not just "group"
        #legendgrouptitle_text="Phaseblocks",
        x=phaseblocks_positions_hp1,
        y=haplotype_1_phaseblocks_values,
        name="HP-1",
        text=phaseblocks_positions_hp1,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'gray', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['dimgray', 'darkgray', 'white']*len(phaseblocks_positions_hp1),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
        #legendgrouptitle_text="Phaseblocks",
    ))

    fig.add_trace(go.Scatter(
        #legendgroup="group3",
        x=phaseblocks_positions_hp2,
        y=haplotype_2_phaseblocks_values,
        name="HP-2",
        text=phaseblocks_positions_hp2,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'green', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['darkgreen', 'limegreen', 'white']*len(phaseblocks_positions_hp2),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
    ))

def add_scatter_trace_phaseblocks(fig, phaseblocks_positions, haplotype_1_phaseblocks_values, haplotype_2_phaseblocks_values):
    fig.add_trace(go.Scatter(
        #legendgroup="group3",  # this can be any string, not just "group"
        #legendgrouptitle_text="Phaseblocks",
        x=phaseblocks_positions,
        y=haplotype_1_phaseblocks_values,
        name="HP-1",
        text=phaseblocks_positions,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'gray', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['dimgray', 'darkgray', 'white']*len(phaseblocks_positions),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
        #legendgrouptitle_text="Phaseblocks",
    ))

    fig.add_trace(go.Scatter(
        #legendgroup="group3",
        x=phaseblocks_positions,
        y=haplotype_2_phaseblocks_values,
        name="HP-2",
        text=phaseblocks_positions,
        #yaxis="y5",
        line = dict(shape = 'spline', color = 'green', width= 1, dash = 'solid'),
        mode='lines+markers',
        marker={"size": 5},
        opacity=0.5,
        marker_color=['darkgreen', 'limegreen', 'white']*len(phaseblocks_positions),
        showlegend=True,
        marker_symbol='diamond-wide',
        hoverinfo = "x+name+y+text",
        #legendgroup="group2",
    ))

def slice_list_sums(input):
    res, last = [[]], None
    for x in sorted(input):
        if last is None or abs(last - x) <= 1.5:
            res[-1].append(x)
        else:
            res.append([x])
        last = x
    first = [res[0][0]]
    print(res)

    return first + [sum(sub_list) / len(sub_list) for sub_list in res[1:]]