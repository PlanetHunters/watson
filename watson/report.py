# import os
import datetime
import logging
import os
import pathlib
from astropy.coordinates import Angle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, PageBreak, \
    Image, Table, TableStyle, ListFlowable
from os import path
from astropy import units as u
import pandas as pd
import numpy as np
from pdf2image import pdf2image


width, height = A4
resources_dir = path.join(path.dirname(__file__))


class Report:
    LOGO_IMAGE = resources_dir + "/resources/images/watson.png"

    def __init__(self, data_dir, file_name, object_id, ra, dec, t0, period, duration, depth, transit_t0s_list,
                 summary_list_t0s_indexes, v, j, h, k, with_tpfs=True, is_summary=False):
        self.data_dir = data_dir
        self.file_name = file_name
        self.object_id = object_id
        self.ra = ra
        self.dec = dec
        self.t0 = t0
        self.period = period
        self.duration = duration
        self.depth = depth
        self.transit_t0s_list = transit_t0s_list
        self.summary_list_t0s_indexes = summary_list_t0s_indexes
        self.v = v
        self.j = j
        self.h = h
        self.k = k
        self.with_tpfs = with_tpfs
        self.is_summary = is_summary

    @staticmethod
    def row_colors(df, table_object):
        data_len = len(df)
        for each in range(1, data_len + 1):
            if each % 2 == 1:
                bg_color = colors.whitesmoke
            else:
                bg_color = colors.lightgrey
            table_object.setStyle(TableStyle([('BACKGROUND', (0, each), (-1, each), bg_color)]))

    @staticmethod
    def metrics_row_colors(df, table_object):
        for index, row in df.iterrows():
            if not isinstance(row['passed'], str):
                if np.isnan(row['passed']):
                    bg_color = colors.yellow
                elif row['passed'] is False or row['passed'] == 0:
                    bg_color = colors.red
                elif row['passed'] is True or row['passed'] == 1:
                    bg_color = colors.lightgreen
                table_index = index + 1
                table_object.setStyle(TableStyle([('BACKGROUND', (0, table_index), (-1, table_index), bg_color),
                                              ('FONTSIZE', (0, table_index), (-1, table_index), 9),
                                              ('TOPPADDING', (0, table_index), (-1, table_index), 1),
                                              ('BOTTOMPADDING', (0, table_index), (-1, table_index), 1)]))

    def create_header(self, canvas, doc):
        canvas.saveState()

        # Logo:
        canvas.drawImage(self.LOGO_IMAGE, x=1.5 * cm, y=26.8 * cm, height=2 * cm, width=2 * cm,
                         preserveAspectRatio=True)

        # Title:
        object_id_text = 'WATSON Transits Validation Report: %s' % self.object_id
        canvas.setFont(psfontname="Helvetica", size=12)
        canvas.drawRightString(x=13.5 * cm, y=27.5 * cm, text=object_id_text)
        if doc.page == 1:
            object_id_text = '%s TRANSITS VALIDATION REPORT' % self.object_id
            canvas.setFont(psfontname="Helvetica-Bold", size=23)
            canvas.drawCentredString(x=10 * cm, y=25.5 * cm, text=object_id_text)

        # Report date:
        report_date = datetime.datetime.now().strftime("%a, %d %B %Y, %H:%M:%S")
        report_date_text = '%s' % report_date

        canvas.setFont("Helvetica", 9)
        canvas.drawRightString(20.5 * cm, 28 * cm, report_date_text)

        canvas.restoreState()

    def create_footer(self, canvas, doc):
        canvas.saveState()

        # if doc.page == 1:
        #     # Footer con superíndice:
        #     textobject = canvas.beginText()
        #     textobject.setTextOrigin(1.8 * cm, 2.1 * cm)
        #     textobject.setFont("Helvetica", 5)
        #     textobject.setRise(5)
        #     textobject.textOut('1 ')
        #     textobject.setRise(0)
        #     textobject.setFont("Helvetica", 7)
        #     pie_pagina = 'Three possible observability values are defined: 1 - Entire transit is required, ' \
        #                  '0.5 - Transit midtime and either ingress or egress at least are required,\n' \
        #                  '0.25 - Only ingress or egress are required, with moon constraints of % sº as minimum ' \
        #                  'distance for new moon and % sº as minimum distance for full moon\n' \
        #                  'and for the observatories listed in the Table 2.' % (self.min_dist, self.max_dist)
        #
        #     for line in pie_pagina.splitlines():
        #         textobject.textLine(line)
        #
        #     canvas.drawText(textobject)

        # Powered by:
        page = "Powered by ReportLab"
        canvas.setFont("Helvetica", 9)
        canvas.drawRightString(7 * cm, 0.5 * cm, page)

        # Page:
        page = "Page %s" % doc.page
        canvas.setFont("Helvetica", 9)
        canvas.drawRightString(20.5 * cm, 0.5 * cm, page)

        canvas.restoreState()

    def create_report(self):
        # Styles to be used
        styles = getSampleStyleSheet()
        styles.add(ParagraphStyle(name="ParagraphAlignCenter", alignment=TA_CENTER))
        styles.add(ParagraphStyle(name="ParagraphAlignJustify", alignment=TA_JUSTIFY))
        styles.wordWrap = 'LTR'
        table_style = TableStyle([('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                                  ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                                  ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
                                  ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
                                  ('FONTSIZE', (0, 0), (-1, -1), 10),
                                  ])
        table_style_small = TableStyle([('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                                        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                                        ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
                                        ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
                                        ('FONTSIZE', (0, 0), (-1, -1), 8),
                                        ])
        # Content:
        story = [Spacer(1, 75)]
        section = 1
        story.append(Paragraph("Section " + str(section) + ": Introduction", styles["Heading1"]))
        introduction = '<font name="HELVETICA" size="9">This document is created by the WATSON report generator (' \
                       '<a href="https://github.com/PlanetHunters/watson" color="blue">https://github.com/PlanetHunters/watson</a>) ' \
                       'and focuses on the target star %s.</font>' % self.object_id
        story.append(Paragraph(introduction, styles["ParagraphAlignJustify"]))

        story.append(Spacer(1, 30))

        # Generamos la tabla 1 con los parámetros:
        tabla1_data = [['RA (deg)', 'Dec (deg)', 'V (mag)', 'J (mag)', 'H (mag)', 'K (mag)'],
                       [Angle(self.ra, u.deg).to_string(unit=u.hourangle, sep=':',
                                                        precision=2) if self.ra is not None else '-',
                        Angle(self.dec, u.deg).to_string(unit=u.deg, sep=':',
                                                         precision=2) if self.dec is not None else '-',
                        round(self.v, 2) if self.v is not None else '-',
                        round(self.j, 2) if self.j is not None else '-',
                        round(self.h, 2) if self.h is not None else '-',
                        round(self.k, 2) if self.k is not None else '-']]
        table1_colwidth = [3.5 * cm, 3.5 * cm, 2 * cm, 2 * cm, 2 * cm, 2 * cm]
        table1_number_rows = len(tabla1_data)
        tabla1 = Table(tabla1_data, table1_colwidth, table1_number_rows * [0.75 * cm])
        tabla1.setStyle(table_style)
        # Le damos el estilo alternando colores de filas:
        Report.row_colors(tabla1_data, tabla1)
        story.append(tabla1)
        table_no = 1
        table1_descripcion = '<font name="HELVETICA" size="9"><strong>Table ' + str(table_no) + ': </strong>\
                        The proposed target parameters.</font>'
        story.append(Spacer(1, 5))
        story.append(Paragraph(table1_descripcion, styles["ParagraphAlignCenter"]))
        story.append(Spacer(1, 15))
        # Generamos la tabla 2 con los parámetros:
        tabla2_data = [['T0 (d)', 'Period (d)', 'Duration (h)', 'Depth (ppt)'],
                       [round(self.t0, 4),
                        round(self.period, 4),
                        round(self.duration / 60, 2),
                        round(self.depth, 3)]]
        table2_colwidth = [4 * cm, 4 * cm, 3.5 * cm, 3.5 * cm]
        table2_number_rows = len(tabla2_data)
        tabla2 = Table(tabla2_data, table2_colwidth, table2_number_rows * [0.75 * cm])
        tabla2.setStyle(table_style)
        Report.row_colors(tabla2_data, tabla2)
        story.append(tabla2)
        story.append(Spacer(1, 5))
        table_no = table_no + 1
        table2_descripcion = '<font name="HELVETICA" size="9"><strong>Table ' + str(table_no) + ': </strong>' \
                             'The candidate parameters.</font>'
        story.append(Paragraph(table2_descripcion, styles["ParagraphAlignCenter"]))
        story.append(Spacer(1, 15))
        metrics_file = self.data_dir + "/metrics.csv"
        iatson_original_predictions_file = f'{self.data_dir}/iatson_predictions.csv'
        iatson_averages_file = f'{self.data_dir}/iatson_averages.csv'
        iatson_branches_file = f'{self.data_dir}/iatson_explain_branches.csv'
        iatson_values_file = f'{self.data_dir}/iatson_explain_values.csv'
        gpt_file = self.data_dir + '/gpt.csv'
        triceratops_validation_file = self.data_dir + '/triceratops/validation.csv'
        table_data = [['Metric', 'Value']]
        explainability_branches_table_data = [['Metric', 'Value']]
        explainability_values_table_data = [['Metric', 'Value']]
        metrics_df = pd.DataFrame(columns=['metric', 'score', 'passed'])
        if os.path.exists(iatson_original_predictions_file):
            averages_df = pd.read_csv(iatson_averages_file)
            branches_df = pd.read_csv(iatson_branches_file)
            values_df = pd.read_csv(iatson_values_file)
            score_average = averages_df.loc[0, 'prediction_value_cal_mean']
            score_std = averages_df.loc[0, 'prediction_value_cal_std']
            # if score_average + score_std > 0.961:
            #     passed = True
            if score_average + score_std >= 0.646:
                passed = True
            elif score_average + score_std <= 0.056:
                passed = False
            else:
                passed = np.nan
            metrics_df = pd.concat([metrics_df, pd.DataFrame.from_dict(
                {"metric": ["WATSON-NET"], 'score': [round(score_average, 4)], 'passed': [passed]}, orient='columns')], ignore_index=True)
            if score_std < 0.015:
                passed = True
            elif score_std < 0.1:
                passed = np.nan
            metrics_df = pd.concat([metrics_df, pd.DataFrame.from_dict(
                {"metric": ["WATSON-NET err"], 'score': [round(score_std, 4)], 'passed': [passed]}, orient='columns')], ignore_index=True)
            for index, row in branches_df.iterrows():
                explainability_branches_table_data.append([row['object_id'], round(row['prediction_value_cal_mean'], 3)])
            for index, row in values_df.iterrows():
                explainability_values_table_data.append([row['object_id'], round(row['prediction_value_cal_mean'], 6)])
        if os.path.exists(gpt_file):
            gpt_df = pd.read_csv(gpt_file)
            score_gpt = gpt_df.loc[0, 'prediction']
            content_gpt = gpt_df.loc[0, 'content']
            if score_gpt == 1:
                passed = True
            elif score_gpt == 0:
                passed = False
            else:
                passed = np.nan
            metrics_df = pd.concat([metrics_df, pd.DataFrame.from_dict(
                {"metric": ["GPT"], 'score': [score_gpt], 'passed': [passed]},
                orient='columns')], ignore_index=True)
        if os.path.exists(triceratops_validation_file):
            validation_df = pd.read_csv(triceratops_validation_file)
            validation_mean_scenario_df = validation_df.loc[validation_df['scenario'] == 'MEAN']
            mean_fpp = validation_mean_scenario_df.reset_index().loc[0, 'FPP']
            if mean_fpp <= 0.015:
                passed = True
            elif mean_fpp <= 0.5:
                passed = np.nan
            else:
                passed = False
            metrics_df = pd.concat([metrics_df, pd.DataFrame.from_dict(
                {"metric": ["Triceratops_FPP"], 'score': [mean_fpp], 'passed': [passed]},
                orient='columns')], ignore_index=True)
            mean_nfpp = validation_mean_scenario_df.reset_index().loc[0, 'NFPP']
            if mean_nfpp <= 0.001:
                passed = True
            elif mean_nfpp <= 0.1:
                passed = np.nan
            else:
                passed = False
            metrics_df = pd.concat([metrics_df, pd.DataFrame.from_dict(
                {"metric": ["Triceratops_NFPP"], 'score': [mean_nfpp], 'passed': [passed]},
                orient='columns')], ignore_index=True)
        story.append(PageBreak())
        section = section + 1
        story.append(Paragraph("Section " + str(section) + ": Metrics summary", styles["Heading1"]))
        if os.path.exists(metrics_file):
            metrics_df = pd.concat([metrics_df, pd.read_csv(metrics_file)], ignore_index=True)
            metrics_df['passed'] = metrics_df['passed'].replace({'True': 1, 'False': 0})
            metrics_df['passed'] = metrics_df['passed'].replace({True: 1, False: 0})
            metrics_df['passed'] = pd.to_numeric(metrics_df['passed'], errors='coerce')
            for index, metric_row in metrics_df.iterrows():
                table_data.append([metric_row['metric'], round(metric_row['score'], 3)])
            table_colwidth = [4 * cm, 4 * cm, 3.5 * cm]
            table_number_rows = len(table_data)
            table = Table(table_data, table_colwidth, table_number_rows * [0.5 * cm])
            table.setStyle(table_style)
            Report.metrics_row_colors(metrics_df, table)
            story.append(table)
            story.append(Spacer(1, 5))
            table_no = table_no + 1
            table_descripcion = '<font name="HELVETICA" size="9"><strong>Table ' + str(table_no) + ': </strong>' \
                                + 'The results of the numerical tests. Green cells mean acceptable values. Yellow means ' \
                                + 'inconclusive. Red represents problematic metrics.</font>'
            story.append(Paragraph(table_descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
        if os.path.exists(iatson_original_predictions_file):
            story.append(PageBreak())
            section = section + 1
            story.append(Paragraph("Section " + str(section) + ": WATSON-Net explainability", styles["Heading1"]))
            table_colwidth = [9 * cm, 4 * cm, 3.5 * cm]
            table_number_rows = len(explainability_branches_table_data)
            table = Table(explainability_branches_table_data, table_colwidth, table_number_rows * [0.5 * cm])
            table.setStyle(table_style)
            story.append(table)
            story.append(Spacer(1, 5))
            table_no = table_no + 1
            table_descripcion = '<font name="HELVETICA" size="9"><strong>Table ' + str(table_no) + ': </strong>' \
                                + ('Impact of each WATSON-Net neural network branch on the final predictions. A positive '
                                   'value means that the scenario helps to classify the signal as a transiting planet.</font>')
            story.append(Paragraph(table_descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            table_colwidth = [9 * cm, 4 * cm, 3.5 * cm]
            table_number_rows = len(explainability_values_table_data)
            table = Table(explainability_values_table_data, table_colwidth, table_number_rows * [0.5 * cm])
            table.setStyle(table_style)
            story.append(table)
            story.append(Spacer(1, 5))
            table_no = table_no + 1
            table_descripcion = '<font name="HELVETICA" size="9"><strong>Table ' + str(table_no) + ': </strong>' \
                                + ('Impact of the variation of the single metrics in the WATSON-Net predictions. A positive '
                                   'value means that the used value improves the prediction in comparison to the original '
                                   'one, making it closer to a transiting planet classification.</font>')
            story.append(Paragraph(table_descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
        if os.path.exists(gpt_file):
            story.append(PageBreak())
            section = section + 1
            story.append(Paragraph("Section " + str(section) + ": WATSON-Net explainability", styles["Heading1"]))
            story.append(Paragraph('<font name="HELVETICA" size="9">GPT has been enabled to analyze this report. Its output is:</font>', styles["ParagraphAlignJustify"]))
            for gpt_explanation_paragraph in content_gpt.split('\n'):
                story.append(Paragraph('<font name="HELVETICA" size="9">' + gpt_explanation_paragraph + '</font>', styles["ParagraphAlignJustify"]))
            story.append(Spacer(1, 30))
        figure = 1
        # TRICERATOPS
        triceratops_base_path = self.data_dir + '/triceratops/'
        if os.path.exists(triceratops_base_path):
            story.append(PageBreak())
            section = section + 1
            story.append(Paragraph("Section " + str(section) + ": TRICERATOPS results", styles["Heading1"]))
        contrast_curve_file = triceratops_base_path + "contrast_curve.png"
        if os.path.exists(contrast_curve_file):
            story.append(Image(contrast_curve_file, width=16 * cm, height=10 * cm))
            descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(
                figure) + ': </strong>Contrast curve for target.</font>'
            story.append(Spacer(1, 5))
            story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            figure = figure + 1
        all_files = os.listdir(triceratops_base_path)
        fov_file = triceratops_base_path + '/fov.png'
        for file in all_files:
            if file.startswith('field_'):
                images = pdf2image.convert_from_path(triceratops_base_path + '/' + file)
                for i in range(len(images)):
                    # Save pages as images in the pdf
                    images[i].save(fov_file, 'PNG')
                break
        story.append(Image(fov_file, width=16 * cm, height=7 * cm))
        descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(
            figure) + ': </strong>Nearby stars for target and its aperture</font>'
        story.append(Spacer(1, 5))
        story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
        story.append(Spacer(1, 15))
        figure = figure + 1
        validation_file = triceratops_base_path + "/validation_scenarios.csv"
        if os.path.exists(validation_file):
            table_data = [['ID', 'scenario', 'M_s', 'R_s', 'P_orb', 'inc', 'b', 'ecc', 'w', 'R_p', 'M_EB', 'R_EB',
                           'prob']]
            metrics_df = pd.read_csv(validation_file)
            for index, metric_row in metrics_df.iterrows():
                table_data.append([str(metric_row['ID']),
                                   metric_row['scenario'],
                                   round(metric_row['M_s'], 2),
                                   round(metric_row['R_s'], 2),
                                   round(metric_row['P_orb'], 2),
                                   round(metric_row['inc'], 2),
                                   round(metric_row['b'], 2),
                                   round(metric_row['ecc'], 2),
                                   round(metric_row['w'], 2),
                                   round(metric_row['R_p'], 2),
                                   round(metric_row['M_EB'], 2),
                                   round(metric_row['R_EB'], 2),
                                   round(metric_row['prob'], 6)])
            table_colwidth = [2.3 * cm, 2 * cm, 1 * cm, 1 * cm, 1 * cm, 1 * cm, 1 * cm, 1 * cm, 1 * cm, 1 * cm,
                              1 * cm, 1 * cm, 2.5 * cm]
            table_number_rows = len(table_data)
            table = Table(table_data, table_colwidth, table_number_rows * [0.5 * cm])
            table.setStyle(table_style_small)
            Report.row_colors(metrics_df, table)
            story.append(table)
            story.append(Spacer(1, 5))
            table_no = table_no + 1
            table_descripcion = '<font name="HELVETICA" size="9"><strong>Table ' + str(table_no) + ': </strong>' \
                                'TRICERATOPS scenarios attributes and probabilities.</font>'
            story.append(Paragraph(table_descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
        validation_file = triceratops_base_path + "/validation.csv"
        if os.path.exists(validation_file):
            table_data = [['Scenario', 'FPP', 'FPP2', 'FPP3+', 'FPP sys', 'FPP2 sys', 'FPP3+ sys', 'NFPP']]
            metrics_df = pd.read_csv(validation_file)
            for index, metric_row in metrics_df.iterrows():
                table_data.append([metric_row['scenario'],
                                   round(metric_row['FPP'], 6),
                                   round(metric_row['FPP2'], 6),
                                   round(metric_row['FPP3+'], 6),
                                   round(metric_row['FPP_sys'], 6),
                                   round(metric_row['FPP2_sys'], 6),
                                   round(metric_row['FPP3+_sys'], 6),
                                   round(metric_row['NFPP'], 6)])
            table_colwidth = [2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm]
            table_number_rows = len(table_data)
            table = Table(table_data, table_colwidth, table_number_rows * [0.5 * cm])
            table.setStyle(table_style)
            Report.row_colors(metrics_df, table)
            story.append(table)
            story.append(Spacer(1, 5))
            table_no = table_no + 1
            table_descripcion = '<font name="HELVETICA" size="9"><strong>Table ' + str(table_no) + ': </strong>' \
                                'TRICERATOPS validation results.</font>'
            story.append(Paragraph(table_descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
        result_map_file = triceratops_base_path + "/triceratops_map.png"
        if os.path.exists(result_map_file):
            story.append(Image(result_map_file, width=9 * cm, height=9 * cm))
            descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(
                figure) + ': </strong>TRICERATOPS validation map for the mean scenario.</font>'
            story.append(Spacer(1, 5))
            story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            figure = figure + 1
        scenario = 0
        scenarios_file = triceratops_base_path + "/scenario_" + str(scenario) + "_fits.pdf"
        scenarios_png_file = triceratops_base_path + "/scenario_" + str(scenario) + "_fits.png"
        while os.path.exists(scenarios_file):
            images = pdf2image.convert_from_path(scenarios_file)
            for i in range(len(images)):
                # Save pages as images in the pdf
                images[i].save(scenarios_png_file, 'PNG')
            story.append(Image(scenarios_png_file, width=16 * cm, height=20 * cm))
            descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(
                figure) + ': </strong>TRICERATOPS best fits for scenario no. ' + str(scenario) + '.</font>'
            story.append(Spacer(1, 5))
            story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            figure = figure + 1
            scenario = scenario + 1
            scenarios_file = triceratops_base_path + "/scenario_" + str(scenario) + "_fits.pdf"
            scenarios_png_file = triceratops_base_path + "/scenario_" + str(scenario) + "_fits.png"

        # OFFSETS
        # neighbours_file_index = 0
        # neighbours_file = self.data_dir + "/star_nb_" + (neighbours_file_index) + ".png"
        # figure = 1
        # while os.path.exists(neighbours_file):
        #     story.append(Image(neighbours_file, width=16 * cm, height=16 * cm))
        #     descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + ': </strong>' \
        #                          'Nearby stars folded plot with the same period than the candidate.</font>'
        #     story.append(Spacer(1, 5))
        #     story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
        #     story.append(Spacer(1, 15))
        #     neighbours_file_index = neighbours_file_index + 1
        #     neighbours_file = self.data_dir + "/star_nb_" + str(neighbours_file_index) + ".png"
        #     figure = figure + 1
        story.append(PageBreak())
        section = section + 1
        story.append(Paragraph("Section " + str(section) + ": Vetting information", styles["Heading1"]))
        transit_depths_file = self.data_dir + "/transit_depths.png"
        if os.path.exists(transit_depths_file):
            story.append(Image(transit_depths_file, width=14 * cm, height=8 * cm))
            descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + ': </strong>' \
                                                                                            'The candidate single-transits depths plot.</font>'
            story.append(Spacer(1, 5))
            story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            figure = figure + 1
        story.append(Image(self.data_dir + "/odd_even_folded_curves.png", width=10 * cm, height=14 * cm))
        descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + ': </strong>' \
                                                                                        'Above, the candidate folded at its found period for the found epoch. ' \
                                                                                        'Middle, the candidate secondary event plot. ' \
                                                                                        'Bottom, the candidate odd and even curves and their binned difference with its best fitted model.</font>'
        story.append(Spacer(1, 5))
        story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
        story.append(Spacer(1, 15))
        figure = figure + 1
        source_offsets_file = self.data_dir + '/source_offsets.png'
        if os.path.exists(source_offsets_file):
            story.append(Image(source_offsets_file, width=14 * cm, height=21 * cm))
            descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + ': </strong>' \
                                                                                            'Above, the computed target offset (red circle) from the differential image offset (cyan dot)' \
                                                                                            ' and the per-pixel BLS SNR offset (green dot). Middle left, the right ascension centroid shift with binning. ' \
                                                                                            'Middle right, the declination centroid shift with binning. ' \
                                                                                            'Bottom left, optical ghost diagnostic curve for core flux.' \
                                                                                            'Bottom right, optical ghost diagnostic curve for halo flux</font>'
            story.append(Spacer(1, 5))
            story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            figure = figure + 1
        if not self.is_summary:
            for file in sorted(list(pathlib.Path(self.data_dir).glob('folded_tpf_*.png'))):
                story.append(Image(str(file), width=14 * cm, height=22 * cm))
                descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + '' \
                                                                                                ': </strong>Above, the TPF and per-pixel BLS SNR best fits. Bottom left, the per-pixel BLS SNR for each' \
                                                                                                ' pixel. Bottom right, the differential images SNR for each pixel. The target position is represented ' \
                                                                                                'by a red star and the TPF independent source offset is represented by a white plus.</font>'
                story.append(Spacer(1, 5))
                story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
                figure = figure + 1
            story.append(PageBreak())
        introduction = '<font name="HELVETICA" size="9">The next pages will contain each of the single-transits ' \
                       'vetting sheets with the next information: </font>'
        story.append(Paragraph(introduction, styles["ParagraphAlignJustify"]))
        story.append(ListFlowable(
            [Paragraph(s, styles["ParagraphAlignJustify"]) for s in [
                "<b>TOP-LEFT</b>: Plot with single transit photometry found in the analyzed curve (with momentum dumps, if any).",
                "<b>TOP-CENTER</b>: Plot with X-axis data drift (X-axis motion vs X-axis centroid offset) around the transit times.",
                "<b>TOP-RIGHT</b>: Plot with Y-axis data drift (Y-axis motion vs Y-axis centroid offset) around the transit times.",
                "<b>CENTER-LEFT</b>: Plot with SAP for used aperture vs SAP for smaller aperture around the transit times.",
                "<b>CENTER-CENTER</b>: Plot with smaller aperture over used aperture on the target.",
                "<b>CENTER-RIGHT</b>: Plot with single transit photometry found in the analyzed curve around the transit times.",
                "<b>BOTTOM</b>: Plot TPF flux measurements for each pixel around the transit times.",
            ]],
            bulletFormat="%s.",  # ←
        ))
        # Pasamos a la siguiente página:
        story.append(PageBreak())
        width = 17 * cm
        height = 24 * cm if self.with_tpfs else 10 * cm
        for index, transit_time in enumerate(self.transit_t0s_list):
            if self.summary_list_t0s_indexes is None or (self.summary_list_t0s_indexes is not None and
                                                         index in self.summary_list_t0s_indexes):
                image_file = self.data_dir + "/single_transit_" + str(index) + "_T0_" + str(transit_time) + ".png"
                if not os.path.exists(image_file):
                    logging.warning(f"The image {image_file} doesn't exist, skipping.")
                    continue
                image = Image(image_file, width=width, height=height)
                story.append(image)
                table3_descripcion = '<font name="HELVETICA" size="9"><strong>Figure %s: </strong>' \
                                     'The single transit no. %s vetting plots</font>' % (str(figure), str(index))
                story.append(Spacer(1, 5))
                story.append(Paragraph(table3_descripcion, styles["ParagraphAlignCenter"]))
                story.append(PageBreak() if self.with_tpfs or figure % 2 == 0 else Spacer(1, 5))
                figure = figure + 1

        # Construimos el documento:
        global_frame = Frame(1.5 * cm, 1.1 * cm, 18 * cm, 25.4 * cm, id='normal', showBoundary=0)
        global_template = PageTemplate(id='UnaColumna', frames=global_frame,
                                       onPage=self.create_header, onPageEnd=self.create_footer)
        doc = BaseDocTemplate(self.data_dir + "/" + self.object_id + "_" + self.file_name, pagesize=A4,
                              rightMargin=40, leftMargin=40,
                              topMargin=95, bottomMargin=15,
                              pageTemplates=global_template)
        doc.build(story)
