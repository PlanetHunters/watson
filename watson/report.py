# import os
import datetime
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
            if np.isnan(row['passed']):
                bg_color = colors.yellow
            elif row['passed'] == False:
                bg_color = colors.red
            elif row['passed'] == True:
                bg_color = colors.lightgreen
            table_object.setStyle(TableStyle([('BACKGROUND', (0, index + 1), (-1, index + 1), bg_color),
                                              ('FONTSIZE', (0, index + 1), (-1, index + 1), 9),
                                              ('TOPPADDING', (0, index + 1), (-1, index + 1), 1),
                                              ('BOTTOMPADDING', (0, index + 1), (-1, index + 1), 1)]))

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
        # Content:
        story = [Spacer(1, 75)]
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
        table1_descripcion = '<font name="HELVETICA" size="9"><strong>Table 1: </strong>\
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
        table2_descripcion = '<font name="HELVETICA" size="9"><strong>Table 2: </strong>' \
                             'The candidate parameters.</font>'
        story.append(Paragraph(table2_descripcion, styles["ParagraphAlignCenter"]))
        story.append(Spacer(1, 15))
        metrics_file = self.data_dir + "/metrics.csv"
        if os.path.exists(metrics_file):
            table_data = [['Metric', 'Value', 'Passed']]
            metrics_df = pd.read_csv(metrics_file)
            for index, metric_row in metrics_df.iterrows():
                table_data.append([metric_row['metric'],
                                   round(metric_row['score'], 3),
                                   metric_row['passed']])
            table_colwidth = [4 * cm, 4 * cm, 3.5 * cm]
            table_number_rows = len(table_data)
            table = Table(table_data, table_colwidth, table_number_rows * [0.5 * cm])
            table.setStyle(table_style)
            Report.metrics_row_colors(metrics_df, table)
            story.append(table)
            story.append(Spacer(1, 5))
            table_descripcion = '<font name="HELVETICA" size="9"><strong>Table 3: </strong>' \
                                'The results of the numerical tests.</font>'
            story.append(Paragraph(table_descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
        cadences_file = self.data_dir + "/folded_cadences.png"
        figure = 1
        if os.path.exists(cadences_file):
            story.append(Image(cadences_file, width=16 * cm, height=20 * cm))
            descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + ': </strong>' \
                                                                                            'Folded curve for all available cadences.</font>'
            story.append(Spacer(1, 5))
            story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            figure = figure + 1
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
        transit_depths_file = self.data_dir + "/transit_depths.png"
        if os.path.exists(transit_depths_file):
            story.append(Image(transit_depths_file, width=14 * cm, height=8 * cm))
            descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + ': </strong>' \
                                                                                            'The candidate single-transits depths plot.</font>'
            story.append(Spacer(1, 5))
            story.append(Paragraph(descripcion, styles["ParagraphAlignCenter"]))
            story.append(Spacer(1, 15))
            figure = figure + 1
        story.append(Image(self.data_dir + "/odd_even_folded_curves.png", width=16 * cm, height=14 * cm))
        descripcion = '<font name="HELVETICA" size="9"><strong>Figure ' + str(figure) + ': </strong>' \
                                                                                        'Above, the candidate folded at its found period for the found epoch and epoch + P/2.' \
                                                                                        'Middle, the candidate folded at its harmonic for the found epoch and epoch + P.' \
                                                                                        'Bottom, the candidate folded at its subharmonic for the found epoch and epoch + P/2, where ' \
                                                                                        'the candidate has been masked.</font>'
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
                image = Image(self.data_dir + "/single_transit_" + str(index) + "_T0_" + str(transit_time) + ".png",
                              width=width, height=height)
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
