# This file is part of BLASTnBRUSH
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the GPL v3.0 License

import matplotlib.patches
import sys
from PyQt6.QtWidgets import (
    QListWidgetItem,
    QApplication, QMainWindow, QWidget, QTabWidget, QLabel, QVBoxLayout, QHBoxLayout,
    QPushButton, QCheckBox, QListWidget, QTextEdit, QLineEdit, QFileDialog, QSpinBox, QDoubleSpinBox,
    QStatusBar, QFrame, QSpacerItem, QSizePolicy, QDialog, QVBoxLayout, QLabel, QPlainTextEdit, QPushButton
)
from PyQt6.QtWidgets import QProgressDialog
from PyQt6.QtGui import QPixmap, QIcon, QBrush, QColor
from PyQt6.QtCore import Qt, QItemSelectionModel
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches


class ChartWidget(QWidget):
    def __init__(self, title):
        super().__init__()
        layout = QVBoxLayout()
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0) 
        self.label = QLabel(title)
        #layout.addWidget(self.label)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.hist_lengths = []
        self.al_lengths = []
        min_seq_len = 0
        max_seq_len = 0
        min_score = 0
        max_score = 0
        min_eval = 0
        max_eval = 0

############## LIBRARIES #######################
from thread_utils import DeleteDuplicatesWorker
from fasta_utils import load_fasta_file, save_fasta_to_db, fetch_all_sequences, fetch_advanced_filtered_sequences
from draw_utils import draw_length_histogram, draw_bitscore_histogram, draw_evalue_histogram, draw_alength_histogram, draw_identities_histogram, draw_positives_histogram
from move_utils import move_checked_items
from blast_utils import load_prev_database, choose_database, choose_database_n, run_blast,run_blastn,run_blastx, save_sequences_from_blast, run_tblastn, run_tblastx, open_sequence_from_xml, show_blast_file, save_blast_output_as
from blast_view import parse_blast_xml, BlastAlignmentViewer
from filter_worker import FilterWorker
from filter_thread import start_parallel_filtering
from gene_loader import GeneLoaderWorker

class MainWindow(QMainWindow):

    ################# STARTUP ###################################
    def __init__(self):
        super().__init__()
        self.threads = []
        self.workers = []
        self.all_results = []
        min_seq_len = 0
        max_seq_len = 0
        min_score = 0
        max_score = 0
        min_eval = 0
        max_eval = 0
        self.setWindowTitle("BLASTnBRUSH")
        self.setWindowIcon(QIcon("Resources/BLASTnBRUSH.png"))

        self.status = QStatusBar()
        self.setStatusBar(self.status)

        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        self.cleaner_tab = QWidget()
        self.blast_tab = QWidget()
        
        self.tabs.addTab(self.blast_tab, "BLAST")
        self.tabs.addTab(self.cleaner_tab, "Cleaner")
        
        self.setup_cleaner_tab()
        self.setup_blast_tab()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< UI >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ########## VERTICAL SEPARATOR #################
    def make_separator(self):
        frame = QFrame()
        frame.setFrameShape(QFrame.Shape.VLine)
        frame.setFrameShadow(QFrame.Shadow.Sunken)
        frame.setFixedWidth(2)
        return frame
    ########## HORIZONTAL SEPARATOR #################
    def make_hor_separator(self):
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFrameShadow(QFrame.Shadow.Sunken)
        separator.setFixedHeight(2)
        separator.setStyleSheet("color: gray; background-color: gray;")
        return separator
    ########## HORIZONTAL SPACE #################
    def make_hor_space(self):
        hor_space = QFrame()
        hor_space.setFrameShape(QFrame.Shape.HLine)
        hor_space.setFrameShadow(QFrame.Shadow.Plain)
        hor_space.setFixedHeight(4)
        return hor_space
           
    ################
    def update_progress_text(self, text):
        self.progress_dialog.setLabelText(text)
    ################
    def update_progress_text(self, text):
        self.progress_dialog.setLabelText(text)
    ################
    def update_progress_value(self, value):
        self.progress_dialog.setValue(value)
        self.progress_dialog.setLabelText(f"Progress: {value}%")
        
    ##################
    def setup_cleaner_tab(self):
        main_layout = QHBoxLayout()

        ############## Left panel ################
        outer_left_layout = QVBoxLayout()
        outer_left_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        inner_left_widget = QWidget()
        left_layout = QVBoxLayout(inner_left_widget)
        #------------------
        self.label_analyse = QLabel("Analyse FASTA")
        font_bold = self.label_analyse.font()
        font_bold.setBold(True)
        font_bold.setPointSize(14)
        self.label_analyse.setFont(font_bold)
        #------------------
        self.label_file = QLabel("File:")
        self.label_file_name = QLabel("None")
        self.load_button = QPushButton("Choose FASTA Files")
        self.load_button.setFixedSize(150, 40)
        self.load_button.clicked.connect(self.load_file)
        #------------------
        self.show_files_button = QPushButton("Show Files")
        self.show_files_button.setFixedSize(150, 30)
        self.show_files_button.clicked.connect(self.show_loaded_files)
        self.show_files_button.setEnabled(False)
        #------------------
        self.label_checkdup = QLabel("Check duplicates in")
        font_bold2 = self.label_analyse.font()
        font_bold2.setBold(True)
        font_bold2.setPointSize(12)
        self.label_checkdup.setFont(font_bold2)
        #------------------
        checkbox_row = QHBoxLayout()
        self.name_box = QCheckBox("Names")
        self.name_box.setChecked(False)
        #------------------
        self.sequence_box = QCheckBox("Sequences")
        self.sequence_box.setChecked(True)
        #------------------        
        self.show_removed_button = QPushButton("Show Removed")
        #------------------    
        self.show_duplicates_button = QPushButton("Show Duplicates")
        self.show_duplicates_button.setFixedSize(150, 30)
        self.show_duplicates_button.clicked.connect(self.show_duplicate_sequences)
        self.show_duplicates_button.setEnabled(False)
        #------------------        
        self.show_removed_button.setFixedSize(150, 30)
        self.show_removed_button.clicked.connect(self.show_removed_sequences)
        self.show_removed_button.setEnabled(False)
        #------------------
        checkbox_row.addWidget(self.name_box)
        checkbox_row.addWidget(self.sequence_box)
        checkbox_row.addWidget(self.show_removed_button)
        #------------------
        self.clean_button = QPushButton("Clean/Analyse")
        self.clean_button.setFixedSize(150, 40)
        self.clean_button.clicked.connect(self.clean_file)
        self.clean_button.setEnabled(False)
        #------------------
        self.label_filter = QLabel("Analyse Parameters")
        self.label_filter.setFont(font_bold2)
        #------------------  
        self.protein_checkbox_obj = QCheckBox("Protein")
        self.protein_checkbox_obj.setChecked(True)
        self.protein_checkbox_obj.setEnabled(False)
        self.protein_checkbox_obj.stateChanged.connect(self.update_atg_checkbox_text)
        #------------------
        self.gene_checkbox_obj = QCheckBox("Gene")
        self.gene_checkbox_obj.setEnabled(False)
        self.gene_checkbox_obj.stateChanged.connect(self.update_atg_checkbox_text)
        #------------------
        self.atg_checkbox = QCheckBox("Check M")
        self.atg_checkbox.setEnabled(False)
        #------------------
        filter_checkbox_row = QHBoxLayout()
        filter_checkbox_row.addWidget(self.protein_checkbox_obj)
        filter_checkbox_row.addWidget(self.gene_checkbox_obj)
        filter_checkbox_row.addWidget(self.atg_checkbox)
        #------------------
        self.length_checkbox_obj = QCheckBox("Length")
        self.length_checkbox_obj.setEnabled(False)
        #------------------
        
        #------------------
        self.len_box_low = QSpinBox()
        self.len_box_low.setMaximum(900000)
        self.len_box_low.setFixedWidth(70)
        self.len_box_hi = QSpinBox()
        self.len_box_hi.setMaximum(900000)
        self.len_box_hi.setFixedWidth(70)
        #------------------
        length_row = QHBoxLayout()
        length_row.addWidget(self.length_checkbox_obj)
        length_row.addStretch()
        length_row.addWidget(QLabel("MIN:"))
        length_row.addWidget(self.len_box_low)
        length_row.addWidget(QLabel("MAX:"))
        length_row.addWidget(self.len_box_hi)
        #self.label_min = QLabel("MIN")
        #self.label_max = QLabel("MAX")      
        #------------------
        self.score_checkbox_obj = QCheckBox("Score")
        self.score_checkbox_obj.setEnabled(False)
        #------------------      
        self.score_box_low = QSpinBox()
        self.score_box_low.setMaximum(900000)
        self.score_box_low.setFixedWidth(70)
        self.score_box_hi = QSpinBox()
        self.score_box_hi.setMaximum(900000)
        self.score_box_hi.setFixedWidth(70)
        #------------------
        score_row = QHBoxLayout()
        score_row.addWidget(self.score_checkbox_obj)
        score_row.addStretch()
        score_row.addWidget(QLabel("MIN:"))
        score_row.addWidget(self.score_box_low)
        score_row.addWidget(QLabel("MAX:"))
        score_row.addWidget(self.score_box_hi)
        #------------------      
        self.eval_checkbox_obj = QCheckBox("E-Value")
        self.eval_checkbox_obj.setEnabled(False)
        #------------------
        
        self.eval_box_low = QSpinBox()
        self.eval_box_low.setMaximum(900000)
        self.eval_box_low.setFixedWidth(70)
        self.eval_box_hi = QSpinBox()
        self.eval_box_hi.setMaximum(900000)
        self.eval_box_hi.setFixedWidth(70)
        #------------------
        eval_row = QHBoxLayout()
        eval_row.addWidget(self.eval_checkbox_obj)
        eval_row.addStretch()
        eval_row.addWidget(QLabel("10^-"))
        eval_row.addWidget(self.eval_box_low)
        eval_row.addWidget(QLabel("10^-"))
        eval_row.addWidget(self.eval_box_hi)
        #------------------       
        self.alength_checkbox_obj = QCheckBox("Al.Length")
        self.alength_checkbox_obj.setEnabled(False) 
        #------------------ 
        self.alength_box_low = QSpinBox()
        self.alength_box_low.setMaximum(900000)
        self.alength_box_low.setFixedWidth(70)
        self.alength_box_hi = QSpinBox()
        self.alength_box_hi.setMaximum(900000)
        self.alength_box_hi.setFixedWidth(70)
        #------------------ 
        alength_row = QHBoxLayout()
        alength_row.addWidget(self.alength_checkbox_obj)
        alength_row.addStretch()
        alength_row.addWidget(QLabel("MIN:"))
        alength_row.addWidget(self.alength_box_low)
        alength_row.addWidget(QLabel("MAX:"))
        alength_row.addWidget(self.alength_box_hi)
        
        #------------------
        self.identities_checkbox_obj = QCheckBox("Identity")
        self.identities_checkbox_obj.setEnabled(False) 
        
        #------------------
        self.identities_box_low = QSpinBox()
        self.identities_box_low.setMaximum(900000)
        self.identities_box_low.setFixedWidth(70)
        self.identities_box_hi = QSpinBox()
        self.identities_box_hi.setMaximum(900000)
        self.identities_box_hi.setFixedWidth(70)
        #------------------
        identities_row = QHBoxLayout()
        identities_row.addWidget(self.identities_checkbox_obj)
        identities_row.addStretch()
        identities_row.addWidget(QLabel("MIN:"))
        identities_row.addWidget(self.identities_box_low)
        identities_row.addWidget(QLabel("MAX:"))
        identities_row.addWidget(self.identities_box_hi)
        
        #------------------        
        self.positives_checkbox_obj = QCheckBox("Similarity")
        self.positives_checkbox_obj.setEnabled(False) 
        #------------------
        self.positives_box_low = QSpinBox()
        self.positives_box_low.setMaximum(900000)
        self.positives_box_low.setFixedWidth(70)
        self.positives_box_hi = QSpinBox()
        self.positives_box_hi.setMaximum(900000)
        self.positives_box_hi.setFixedWidth(70)
        #------------------
        positives_row = QHBoxLayout()
        positives_row.addWidget(self.positives_checkbox_obj)
        positives_row.addStretch()
        positives_row.addWidget(QLabel("MIN:"))
        positives_row.addWidget(self.positives_box_low)
        positives_row.addWidget(QLabel("MAX:"))
        positives_row.addWidget(self.positives_box_hi)     
        #------------------
        
        
        #------------------
        self.label_status = QLabel("Status:")
        #------------------
        self.label_status_text = QLabel("No file selected")
        #------------------
        self.label_removed = QLabel("Removed 0 sequences")
        #------------------
        self.label_inout = QLabel("Input: 0 Output: 0")
        #------------------
        self.amount_label = QLabel("Records: 0")
        #------------------
        self.amount2_label = QLabel("Records: 0")
        #------------------
        self.label_copyright = QLabel("© Liszka, Marcisz, Stołowski")
        #------------------
        self.analyze_button = QPushButton("Filter")
        self.analyze_button.clicked.connect(self.filter_all)
        self.analyze_button.setFixedSize(150, 40)
        self.analyze_button.setEnabled(False)
        #------------------
        self.reset_button = QPushButton("RESET")
        self.reset_button.setFixedSize(150, 40)
        self.reset_button.clicked.connect(self.reset_filters)
        self.reset_button.setEnabled(False)
        #------------------
        self.save_fasta_button = QPushButton("Save FASTA")
        self.save_fasta_button.clicked.connect(self.save_selected_to_fasta)
        self.save_fasta_button.setFixedSize(150, 40)
        self.save_fasta_button.setEnabled(False)
        #------------------
        #------------------
        #------------------
        left_layout.addWidget(self.label_analyse, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.load_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.show_files_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        #left_layout.addWidget(self.label_file_name, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.make_hor_separator())
        left_layout.addWidget(self.label_checkdup, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addLayout(checkbox_row)
        left_layout.addWidget(self.clean_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_status_text, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_removed, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_inout, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.make_hor_separator())
        left_layout.addWidget(self.label_filter, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addLayout(filter_checkbox_row)
        left_layout.addLayout(length_row)  
        left_layout.addLayout(score_row)  
        left_layout.addLayout(eval_row) 
        left_layout.addLayout(alength_row)
        left_layout.addLayout(identities_row)
        left_layout.addLayout(positives_row)
        left_layout.addWidget(self.make_hor_space())  
        left_layout.addWidget(self.analyze_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.reset_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.save_fasta_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.make_hor_separator())
        #------------------
        self.image_label = QLabel()
        self.image_label.setPixmap(QPixmap("resources/BLASTnBRUSH.png").scaledToWidth(30, Qt.TransformationMode.SmoothTransformation))
        #left_layout.addWidget(self.image_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_copyright, alignment=Qt.AlignmentFlag.AlignHCenter)
        #------------------
        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))
        #------------------     
        outer_left_layout.addWidget(inner_left_widget)
        #------------------
        main_layout.addLayout(outer_left_layout, 1)
        main_layout.addWidget(self.make_separator())
        main_layout.addLayout(self.setup_middle_panel(), 5)
        main_layout.addWidget(self.make_separator())
        main_layout.addLayout(self.setup_right_panel(), 4)
        #------------------
        self.cleaner_tab.setLayout(main_layout)
        
        
     ########## MIDDLE PANEL #################
    def setup_middle_panel(self):
        layout = QVBoxLayout()
        #------------------
        self.label_genes = QLabel("Sequences")
        font_bold = self.label_genes.font()
        font_bold.setBold(True)
        font_bold.setPointSize(14)   
        self.label_genes.setFont(font_bold)
        #------------------
        self.genes_list = QListWidget()
        self.genes_list.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        self.genes_list.setStyleSheet("""
        QListWidget::item:selected {
        background: CornflowerBlue;
        color: white;
        }""")
        self.genes_list.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.genes_list.customContextMenuRequested.connect(self.show_sequence_context_menu)
        #------------------
        self.select_all_checkbox = QCheckBox("Select all")
        self.select_all_checkbox.stateChanged.connect(self.toggle_select_all_genes)
        self.label_selected_names = QLabel("Selected names")
        #------------------
        self.selected_names_listbox = QListWidget()
        self.selected_names_listbox.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        self.selected_names_listbox.setStyleSheet("""
        QListWidget::item:selected {
        background: CornflowerBlue;
        color: white;
        }""")
        #------------------
        self.filter_select_all_checkbox = QCheckBox("Select all")
        self.filter_select_all_checkbox.stateChanged.connect(self.toggle_select_all_selected)
        #------------------
        self.identity_layout = QHBoxLayout()
        self.identity_layout.setSpacing(0)
        self.identity_layout.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.label_similarity = QLabel("Identity [%]")
        #------------------
        self.accuracy_box = QSpinBox()
        self.accuracy_box.setFixedWidth(70)
        self.accuracy_box.setMaximum(100)
        self.accuracy_box.setValue(100)
        #------------------
        self.identity_layout.addWidget(self.label_similarity)
        self.identity_layout.addWidget(self.accuracy_box)
        #------------------
        self.filter_textbox = QTextEdit()
        self.filter_textbox.setPlaceholderText("Paste one or more names (Names separate by [ENTER], parts of single name separate by [ , ])")
        self.filter_textbox.setFixedHeight(60)  # lub inna wysokość
        #------------------
        self.sequence_textbox = QLineEdit()
        #------------------
        self.add_button = QPushButton("ADD")
        self.add_button.clicked.connect(self.add_selected_genes)
        #------------------
        self.remove_button = QPushButton("Remove")
        self.remove_button.clicked.connect(self.remove_selected_items)
        self.reset_button.setFixedSize(150, 40)
        self.reset_button.clicked.connect(self.reset_filters)
        #------------------
        
        self.info1_button = QPushButton("?")
        self.info1_button.clicked.connect(self.info1)
        self.info1_button.setFixedSize(30, 20)
        #------------------
        
        layout.addWidget(self.label_genes, alignment=Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.genes_list)
        #------------------
        sel_all_box = QHBoxLayout()
        sel_all_box.addWidget(self.amount_label)
        sel_all_box.addStretch()
        sel_all_box.addWidget(self.select_all_checkbox)
        layout.addLayout(sel_all_box)
        #------------------
        info1_box = QHBoxLayout()
        info1_box.addWidget(QLabel("Filter by names:"))
        info1_box.addWidget(self.info1_button)
        layout.addLayout(info1_box)
        layout.addWidget(self.filter_textbox)
        layout.addWidget(QLabel("Filter by sequence:"))
        layout.addWidget(self.sequence_textbox)
        layout.addLayout(self.identity_layout)
        layout.addWidget(self.add_button)
        layout.addWidget(self.make_hor_separator())
        layout.addWidget(self.label_selected_names)
        layout.addWidget(self.selected_names_listbox)
        layout.addWidget(self.amount2_label)
        layout.addWidget(self.filter_select_all_checkbox)
        layout.addWidget(self.remove_button)
        
        return layout
        #------------------
        #------------------
    ############## RIGHT PANEL ##########################
    def setup_right_panel(self):

        self.label_histograms = QLabel("Histograms")
        font_bold = self.label_histograms.font()
        font_bold.setBold(True)
        font_bold.setPointSize(14)
        self.label_histograms.setFont(font_bold)
        #------------------
        layout = QVBoxLayout()
        histograms_buttons = QHBoxLayout()
        histograms_buttons.setAlignment(Qt.AlignmentFlag.AlignCenter)
        histograms_buttons2 = QHBoxLayout()
        histograms_buttons2.setAlignment(Qt.AlignmentFlag.AlignCenter)
        histograms_buttons3 = QHBoxLayout()
        histograms_buttons3.setAlignment(Qt.AlignmentFlag.AlignCenter)
        #------------------LENGTH HISTOGRAM BUTTONS------------------------
        self.Length_histogram_button = QPushButton("Length")
        self.Length_histogram_button.setFixedSize(100, 30)
        self.Length_histogram_button.clicked.connect(lambda: draw_length_histogram(self.chart3, self.hist_lengths))
        self.Length_histogram_button2 = QPushButton("Length")
        self.Length_histogram_button2.setFixedSize(100, 30)
        self.Length_histogram_button2.clicked.connect(lambda: draw_length_histogram(self.chart2, self.hist_lengths))
        self.Length_histogram_button3 = QPushButton("Length")
        self.Length_histogram_button3.setFixedSize(100, 30)
        self.Length_histogram_button3.clicked.connect(lambda: draw_length_histogram(self.chart1, self.hist_lengths))
        #------------------SCORE HISTOGRAM BUTTONS------------------------
        self.Score_histogram_button = QPushButton("Score")
        self.Score_histogram_button.setFixedSize(100, 30)
        self.Score_histogram_button.clicked.connect(lambda: draw_bitscore_histogram(self.chart3, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Score_histogram_button2 = QPushButton("Score")
        self.Score_histogram_button2.setFixedSize(100, 30)
        self.Score_histogram_button2.clicked.connect(lambda: draw_bitscore_histogram(self.chart2, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Score_histogram_button3 = QPushButton("Score")
        self.Score_histogram_button3.setFixedSize(100, 30)
        self.Score_histogram_button3.clicked.connect(lambda: draw_bitscore_histogram(self.chart1, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        #------------------EVALUE HISTOGRAM BUTTONS------------------------
        self.Eval_histogram_button = QPushButton("E-value")
        self.Eval_histogram_button.setFixedSize(100, 30)
        self.Eval_histogram_button.clicked.connect(lambda: draw_evalue_histogram(self.chart3, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Eval_histogram_button2 = QPushButton("E-value")
        self.Eval_histogram_button2.setFixedSize(100, 30)
        self.Eval_histogram_button2.clicked.connect(lambda: draw_evalue_histogram(self.chart2, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Eval_histogram_button3 = QPushButton("E-value")
        self.Eval_histogram_button3.setFixedSize(100, 30)
        self.Eval_histogram_button3.clicked.connect(lambda: draw_evalue_histogram(self.chart1, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        #------------------ALIGNMENT HISTOGRAM BUTTONS------------------------
        self.ALength_histogram_button = QPushButton("Al. Length")
        self.ALength_histogram_button.setFixedSize(100, 30)
        self.ALength_histogram_button.clicked.connect(lambda: draw_alength_histogram(self.chart3, self.al_lengths))
        self.ALength_histogram_button2 = QPushButton("Al. Length")
        self.ALength_histogram_button2.setFixedSize(100, 30)
        self.ALength_histogram_button2.clicked.connect(lambda: draw_alength_histogram(self.chart2, self.al_lengths))
        self.ALength_histogram_button3 = QPushButton("Al. Length")
        self.ALength_histogram_button3.setFixedSize(100, 30)
        self.ALength_histogram_button3.clicked.connect(lambda: draw_alength_histogram(self.chart1, self.al_lengths))
        #------------------IDENTITY HISTOGRAM BUTTONS------------------------
        self.Identities_histogram_button = QPushButton("Identity")
        self.Identities_histogram_button.setFixedSize(100, 30)
        self.Identities_histogram_button.clicked.connect(lambda: draw_identities_histogram(self.chart3, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Identities_histogram_button2 = QPushButton("Identity")
        self.Identities_histogram_button2.setFixedSize(100, 30)
        self.Identities_histogram_button2.clicked.connect(lambda: draw_identities_histogram(self.chart2, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Identities_histogram_button3 = QPushButton("Identity")
        self.Identities_histogram_button3.setFixedSize(100, 30)
        self.Identities_histogram_button3.clicked.connect(lambda: draw_identities_histogram(self.chart1, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        #------------------SIMILARITY HISTOGRAM BUTTONS------------------------
        self.Positives_histogram_button = QPushButton("Similarity")
        self.Positives_histogram_button.setFixedSize(100, 30)
        self.Positives_histogram_button.clicked.connect(lambda: draw_positives_histogram(self.chart3, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Positives_histogram_button2 = QPushButton("Similarity")
        self.Positives_histogram_button2.setFixedSize(100, 30)
        self.Positives_histogram_button2.clicked.connect(lambda: draw_positives_histogram(self.chart2, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        self.Positives_histogram_button3 = QPushButton("Similarity")
        self.Positives_histogram_button3.setFixedSize(100, 30)
        self.Positives_histogram_button3.clicked.connect(lambda: draw_positives_histogram(self.chart1, [self.genes_list.item(i).text() for i in range(self.genes_list.count())]))
        #------------------
        self.Length_histogram_button.setEnabled(False)
        self.Score_histogram_button.setEnabled(False)
        self.Eval_histogram_button.setEnabled(False)
        self.ALength_histogram_button.setEnabled(False)
        self.Identities_histogram_button.setEnabled(False)
        self.Positives_histogram_button.setEnabled(False)
        #------------------
        self.Length_histogram_button2.setEnabled(False)
        self.Score_histogram_button2.setEnabled(False)
        self.Eval_histogram_button2.setEnabled(False)
        self.ALength_histogram_button2.setEnabled(False)
        self.Identities_histogram_button2.setEnabled(False)
        self.Positives_histogram_button2.setEnabled(False)
        #------------------
        self.Length_histogram_button3.setEnabled(False)
        self.Score_histogram_button3.setEnabled(False)
        self.Eval_histogram_button3.setEnabled(False)
        self.ALength_histogram_button3.setEnabled(False)
        self.Identities_histogram_button3.setEnabled(False)
        self.Positives_histogram_button3.setEnabled(False)
        #------------------
        histograms_buttons.addWidget(self.Length_histogram_button)
        histograms_buttons.addWidget(self.ALength_histogram_button)
        histograms_buttons.addWidget(self.Score_histogram_button)
        histograms_buttons.addWidget(self.Eval_histogram_button)       
        histograms_buttons.addWidget(self.Identities_histogram_button)
        histograms_buttons.addWidget(self.Positives_histogram_button)
        #------------------
        histograms_buttons2.addWidget(self.Length_histogram_button2)
        histograms_buttons2.addWidget(self.ALength_histogram_button2)
        histograms_buttons2.addWidget(self.Score_histogram_button2)
        histograms_buttons2.addWidget(self.Eval_histogram_button2)       
        histograms_buttons2.addWidget(self.Identities_histogram_button2)
        histograms_buttons2.addWidget(self.Positives_histogram_button2)
        #------------------
        histograms_buttons3.addWidget(self.Length_histogram_button3)
        histograms_buttons3.addWidget(self.ALength_histogram_button3)
        histograms_buttons3.addWidget(self.Score_histogram_button3)
        histograms_buttons3.addWidget(self.Eval_histogram_button3)       
        histograms_buttons3.addWidget(self.Identities_histogram_button3)
        histograms_buttons3.addWidget(self.Positives_histogram_button3)
        #------------------
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.chart1 = ChartWidget("")
        self.chart2 = ChartWidget("")
        self.chart3 = ChartWidget("")
        self.chart4 = ChartWidget("")
        
        layout.addWidget(self.label_histograms, alignment=Qt.AlignmentFlag.AlignHCenter)
        layout.addLayout(histograms_buttons)
        layout.addWidget(self.chart3)
        layout.addLayout(histograms_buttons2)
        layout.addWidget(self.chart2)
        layout.addLayout(histograms_buttons3)
        layout.addWidget(self.chart1)
        return layout
        
        
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
    ######################### LOAD FASTA FILE ########################################
    def load_file(self):
        from fasta_utils import load_fasta_file, save_fasta_to_db
        from PyQt6.QtWidgets import QFileDialog, QMessageBox
        import sqlite3
        import os
        import time

        file_names, _ = QFileDialog.getOpenFileNames(
            self, "Open FASTA Files", "", "FASTA Files (*.fasta *.fa *.txt);;All Files (*)"
        )

        if not file_names:
            return
        
        start = time.time()
        
        self.loaded_files = file_names
        self.label_status_text.setText("Files selected")
        self.clean_button.setEnabled(True)
        self.save_fasta_button.setEnabled(True)
        self.analyze_button.setEnabled(False)
        self.show_files_button.setEnabled(True)

        try:
            with sqlite3.connect("sequences.db") as conn:
                c = conn.cursor()
                c.execute("DROP TABLE IF EXISTS sequences")
                c.execute("CREATE TABLE sequences (id INTEGER PRIMARY KEY AUTOINCREMENT, header TEXT, sequence TEXT)")
                conn.commit()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset database: {e}")
            return

        total_files = len(file_names)
        self.progress_dialog = QProgressDialog("Loading sequences...", "Cancel", 0, total_files, self)
        self.progress_dialog.setWindowTitle("Loading FASTA")
        self.progress_dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress_dialog.setMinimumDuration(0)
        self.progress_dialog.setValue(0)
        self.progress_dialog.show()
        QApplication.processEvents()
               
        for i, file_path in enumerate(file_names):
            save_fasta_to_db(file_path, db_path="sequences.db", append=True)
            self.progress_dialog.setValue(i + 1)
            QApplication.processEvents()

            if self.progress_dialog.wasCanceled():
                break

        self.progress_dialog.close()
        QMessageBox.information(self, "Done", "All FASTA files loaded.")

        self.genes_list.clear()
        self.selected_names_listbox.clear()
        self.min_seq_len = float("inf")
        self.max_seq_len = 0
        
        end = time.time()
        print(f"Load finished in {end - start:.2f} seconds")
        
        
    def show_loaded_files(self):
        from PyQt6.QtWidgets import QMessageBox
        if hasattr(self, "loaded_files") and self.loaded_files:
            msg = QMessageBox(self)
            msg.setWindowTitle("Loaded FASTA Files")
            msg.setText("Files:\n\n" + "\n".join(self.loaded_files))
            msg.exec()
        else:
            QMessageBox.information(self, "Info", "No files loaded.")
    
    def start_loading_genes(self, total_before, total_after):
        self.progress_dialog.setLabelText("Loading genes...")
        self.genes_list.clear()

        self.worker = GeneLoaderWorker()
        self.worker.finished.connect(lambda headers, lengths, al_lengths, stats:
            self.on_genes_loaded(headers, lengths, al_lengths, stats, total_before, total_after))
        self.worker.start()    
            
    ############### RUN CLEANER IN OTHER THREAD #################
    def clean_file(self):
        from PyQt6.QtWidgets import QProgressDialog
        from PyQt6.QtCore import Qt
        import time

        check_name = self.name_box.isChecked()
        check_sequence = self.sequence_box.isChecked()

        self.progress_dialog = QProgressDialog("Starting...", "Cancel", 0, 0, self)
        self.progress_dialog.setWindowTitle("Cleaning Database")
        self.progress_dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress_dialog.setMinimumDuration(0)
        self.progress_dialog.show()

        self.worker = DeleteDuplicatesWorker(check_name, check_sequence)
        self.worker.progress_percent.connect(self.update_progress_value)
        self.worker.progress_text.connect(self.update_progress_text)
        self.worker.finished.connect(self.on_cleaning_finished)
        self.worker.start()
       

    ########################## CLEANING END ############################
    def on_cleaning_finished(self, _unused, total_before, total_after):
        import time
        self.start_time = time.time()
        self.start_loading_genes(total_before, total_after)
    

    def on_genes_loaded(self, headers, lengths, al_lengths, stats, total_before, total_after):
        import time
        import math

        self.genes_list.setUpdatesEnabled(False)
        self.genes_list.addItems(headers)
        self.genes_list.setUpdatesEnabled(True)

        self.hist_lengths = lengths
        self.al_lengths = al_lengths
        self.min_seq_len = stats["min_len"]
        self.max_seq_len = stats["max_len"]

        self.label_removed.setText(f"Removed {total_before - total_after} sequences")
        self.label_inout.setText(f"Input: {total_before} Output: {total_after}")
        self.amount_label.setText(f"Records: {total_after}")
        self.len_box_low.setValue(self.min_seq_len)
        self.len_box_hi.setValue(self.max_seq_len)

        # Enable filters
        self.analyze_button.setEnabled(True)
        self.length_checkbox_obj.setEnabled(True)
        self.protein_checkbox_obj.setEnabled(True)
        self.atg_checkbox.setEnabled(True)
        self.gene_checkbox_obj.setEnabled(True)
        self.len_box_hi.setEnabled(True)
        self.len_box_low.setEnabled(True)
        self.save_fasta_button.setEnabled(True)
        self.clean_button.setEnabled(False)
        self.reset_button.setEnabled(True)
        self.show_removed_button.setEnabled(True)
        self.show_duplicates_button.setEnabled(True)

        # Histogram buttons
        self.Length_histogram_button.setEnabled(True)
        self.Length_histogram_button2.setEnabled(True)
        self.Length_histogram_button3.setEnabled(True)

        if stats["scores"]:
            self.score_checkbox_obj.setEnabled(True)
            self.Score_histogram_button.setEnabled(True)
            self.Score_histogram_button2.setEnabled(True)
            self.Score_histogram_button3.setEnabled(True)
            self.score_box_low.setValue(min(stats["scores"]))
            self.score_box_hi.setValue(max(stats["scores"]))
        else:
            self.score_checkbox_obj.setEnabled(False)

        if stats["e_values"]:
            self.eval_checkbox_obj.setEnabled(True)
            self.Eval_histogram_button.setEnabled(True)
            self.Eval_histogram_button2.setEnabled(True)
            self.Eval_histogram_button3.setEnabled(True)

            # przelicz e-values na potęgi (-log10)
            TINY = 1e-300
            neglogs = [-math.log10(max(e, TINY)) for e in stats["e_values"]]

            # ustaw zakres w potęgach
            self.eval_box_low.setValue(int(min(neglogs)))
            self.eval_box_hi.setValue(int(max(neglogs)))

        else:
            self.eval_checkbox_obj.setEnabled(False)

        if al_lengths:
            self.alength_checkbox_obj.setEnabled(True)
            self.ALength_histogram_button.setEnabled(True)
            self.ALength_histogram_button2.setEnabled(True)
            self.ALength_histogram_button3.setEnabled(True)
            self.alength_box_low.setValue(min(al_lengths))
            self.alength_box_hi.setValue(max(al_lengths))
        else:
            self.alength_checkbox_obj.setEnabled(False)

        if stats["identities"]:
            self.identities_checkbox_obj.setEnabled(True)
            self.Identities_histogram_button.setEnabled(True)
            self.Identities_histogram_button2.setEnabled(True)
            self.Identities_histogram_button3.setEnabled(True)
            self.identities_box_low.setValue(min(stats["identities"]))
            self.identities_box_hi.setValue(max(stats["identities"]))
        else:
            self.identities_checkbox_obj.setEnabled(False)

        if stats["positives"]:
            self.positives_checkbox_obj.setEnabled(True)
            self.Positives_histogram_button.setEnabled(True)
            self.Positives_histogram_button2.setEnabled(True)
            self.Positives_histogram_button3.setEnabled(True)
            self.positives_box_low.setValue(min(stats["positives"]))
            self.positives_box_hi.setValue(max(stats["positives"]))
        else:
            self.positives_checkbox_obj.setEnabled(False)

        # Rysuj histogramy
        #draw_length_histogram(self.chart3, self.hist_lengths)
        #draw_bitscore_histogram(self.chart2, headers)
        #draw_evalue_histogram(self.chart1, headers)

        self.progress_dialog.close()
        print(f"Gene list loaded in {time.time() - self.start_time:.2f} s")

    
    def show_removed_sequences(self):
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QListWidget, QPushButton
        import sqlite3

        conn = sqlite3.connect("duplicates.db")
        cursor = conn.cursor()
        cursor.execute("SELECT header, sequence FROM sequences")
        records = cursor.fetchall()
        conn.close()

        #if not records:
            #return

        dialog = QDialog(self)
        dialog.setWindowTitle("Removed Sequences")
        layout = QVBoxLayout(dialog)

        list_widget = QListWidget()
        for header, seq in records:
            list_widget.addItem(f"{header} [{len(seq)}]")

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)

        layout.addWidget(list_widget)
        layout.addWidget(close_btn)

        dialog.setLayout(layout)
        dialog.resize(700, 500)
        dialog.exec()                                                                                  

    def show_duplicate_sequences(self):
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QListWidget, QPushButton
        import sqlite3

        conn = sqlite3.connect("duplicates.db")
        cursor = conn.cursor()
        cursor.execute("SELECT header, sequence FROM sequences")
        records = cursor.fetchall()
        conn.close()

        dialog = QDialog(self)
        dialog.setWindowTitle("Duplicate Sequences")
        layout = QVBoxLayout(dialog)

        list_widget = QListWidget()
        for header, seq in records:
            list_widget.addItem(f"{header} [{len(seq)}]")

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)

        layout.addWidget(list_widget)
        layout.addWidget(close_btn)

        dialog.setLayout(layout)
        dialog.resize(700, 500)
        dialog.exec()       

        
    ############# ADD GENES TO LIST ############################
    def add_selected_genes(self):
        for i in range(self.genes_list.count()):
            item = self.genes_list.item(i)
            if item.isSelected():
                text = item.text()
                if not any(self.selected_names_listbox.item(j).text() == text for j in range(self.selected_names_listbox.count())):
                    self.selected_names_listbox.addItem(text)
        self.amount2_label.setText(f"Records: {self.selected_names_listbox.count()}")
    
    def prevent_uncheck_name(self, state):
        self.name_box.setChecked(True)    
            
            
    def prevent_uncheck_seq(self, state):        
        self.sequence_box.setChecked(True)
        
    ############# SELECT ALL GENES CHECKBOX ############################
    def toggle_select_all_genes(self, state):
        state == Qt.CheckState.Checked
        for i in range(self.genes_list.count()):
            item = self.genes_list.item(i)
            item.setSelected(state)
    ############# SELECT ALL SELECTED CHECKBOX #########################
    def toggle_select_all_selected(self, state):
        state == Qt.CheckState.Checked
        for i in range(self.selected_names_listbox.count()):
            item = self.selected_names_listbox.item(i)
            item.setSelected(state)
    ############## REMOVE SELECTED FROM LIST ############################
    def remove_selected_items(self):
        for i in reversed(range(self.selected_names_listbox.count())):
            item = self.selected_names_listbox.item(i)
            if item.isSelected():
                self.selected_names_listbox.takeItem(i)
        self.amount2_label.setText(f"Records: {self.selected_names_listbox.count()}")
     
    ################## RIGHT CLICK ON HEADER ############################
    def show_sequence_context_menu(self, position):
        from PyQt6.QtWidgets import QMenu, QMessageBox
        from fasta_utils import fetch_all_sequences

        item = self.genes_list.itemAt(position)
        if item is None:
            return

        menu = QMenu()
        show_action = menu.addAction("Show Sequence")
        action = menu.exec(self.genes_list.viewport().mapToGlobal(position))

        if action == show_action:

            full_text = item.text().split(" [")[0].strip()
            header_text = full_text.split()[0]  

            for header, seq in fetch_all_sequences("cleaned.db"):
                header_main = header.strip().split()[0] 
                if header_main == header_text:
                    self.show_sequence_dialog(header, seq)
                    break
            else:
                QMessageBox.warning(self, "Not found", "Sequence not found in database.")

    ################## SHOW SEQUENCE ############################
    from PyQt6.QtWidgets import QDialog, QVBoxLayout, QLabel, QPlainTextEdit, QPushButton

    def show_sequence_dialog(self, header, sequence):
        dialog = QDialog(self)
        dialog.setWindowTitle("Selected Sequence")

        fasta_text = f"{header}\n{sequence}"
        text_edit = QPlainTextEdit()
        text_edit.setPlainText(fasta_text)
        text_edit.setReadOnly(True)
        
        layout = QVBoxLayout()
        layout.addWidget(text_edit)

        button_layout = QHBoxLayout()

        close_button = QPushButton("Close")
        close_button.clicked.connect(dialog.accept)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

        dialog.setLayout(layout)
        dialog.resize(700, 400)
        dialog.exec()
    
    ####################### FILTERING ############################

    from filter_thread import start_parallel_filtering

    def filter_all(self):
        import time 
        
        start = time.time()
        name_terms = self.filter_textbox.toPlainText().splitlines()
        seq_terms = self.sequence_textbox.text().splitlines()

        filter_args = {
            "name_terms": name_terms,
            "seq_terms": seq_terms,
            "similarity_threshold": self.accuracy_box.value(),

            "check_length": self.length_checkbox_obj.isChecked(),
            "min_len": self.len_box_low.value(),
            "max_len": self.len_box_hi.value(),

            "check_atg": self.atg_checkbox.isChecked(),
            "check_m": self.protein_checkbox_obj.isChecked(),
            "check_atg_dna": self.gene_checkbox_obj.isChecked(),

            "check_score": self.score_checkbox_obj.isChecked(),
            "min_score": self.score_box_low.value(),
            "max_score": self.score_box_hi.value(),

            "check_eval": self.eval_checkbox_obj.isChecked(),
            "min_eval": self.eval_box_low.value(),
            "max_eval": self.eval_box_hi.value(),

            "check_alength": self.alength_checkbox_obj.isChecked(),
            "min_alength": self.alength_box_low.value(),
            "max_alength": self.alength_box_hi.value(),

            "check_identities": self.identities_checkbox_obj.isChecked(),
            "min_identities": self.identities_box_low.value(),
            "max_identities": self.identities_box_hi.value(),

            "check_positives": self.positives_checkbox_obj.isChecked(),
            "min_positives": self.positives_box_low.value(),
            "max_positives": self.positives_box_hi.value(),
        }

        
        self.genes_list.clear()
        self.amount_label.setText("Filtering...")
        self.progress_dialog = QProgressDialog("Filtering...", "Abort", 0, 0, self)
        self.progress_dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress_dialog.show()

        def on_done(results):
            self.progress_dialog.close()
            self.genes_list.clear()
            for header, sequence in results:
                self.genes_list.addItem(f"{header} [{len(sequence)}]")
            self.amount_label.setText(f"Found: {len(results)} records")
            end = time.time()
            print(f"Filtered in {end - start:.2f} seconds")

        start_parallel_filtering(self, filter_args, on_done)
    
    
        
    ################# SAVE FASTA FILE ############################
    def save_selected_to_fasta(self):
        import sqlite3
        from PyQt6.QtWidgets import QFileDialog, QMessageBox, QProgressDialog, QApplication
        from PyQt6.QtCore import Qt, QTimer

        file_name, _ = QFileDialog.getSaveFileName(
            self, "Save Selected FASTA", "", "FASTA Files (*.fasta *.fa);;All Files (*)"
        )

        if not file_name:
            return

        selected_ids = {
            self.selected_names_listbox.item(i).text().split()[0].strip()
            for i in range(self.selected_names_listbox.count())
        }
        
        self.s_progress = QProgressDialog("Saving to FASTA...", "Cancel", 0, 0, self)
        self.s_progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.s_progress.setMinimumDuration(0)
        self.s_progress.setValue(0)
        self.s_progress.show()

        def do_export():
            try:
                with sqlite3.connect("cleaned.db") as conn:
                    c = conn.cursor()
                    c.execute("SELECT COUNT(*) FROM sequences")
                    total_records = c.fetchone()[0]
                    self.s_progress.setMaximum(total_records)

                    matched_count = 0
                    c.execute("SELECT header, sequence FROM sequences")

                    with open(file_name, "w", encoding="utf-8") as f:
                        for i, (header, sequence) in enumerate(c):
                            if self.s_progress.wasCanceled():
                                break

                            header_stripped = header.strip()

                            # obcięcie nazwy do " ["
                            cut_index = header_stripped.find(" [")
                            if cut_index != -1:
                                header_clean = header_stripped[:cut_index]
                            else:
                                header_clean = header_stripped

                            header_id = header_clean.split()[0]

                            if header_id in selected_ids:
                                f.write(f">{header_clean}\n")
                                for j in range(0, len(sequence), 60):
                                    f.write(sequence[j:j+60] + "\n")
                                matched_count += 1

                            if i % 50 == 0 or i == total_records - 1:
                                self.s_progress.setValue(i + 1)
                                QApplication.processEvents()

                self.s_progress.close()

                if matched_count == 0:
                    QMessageBox.warning(self, "Warning", "No matching sequences were saved. Check header formatting.")
                else:
                    QMessageBox.information(self, "Success", f"{matched_count} sequences saved to file!")

            except Exception as e:
                self.s_progress.close()
                QMessageBox.critical(self, "Error", f"Failed to save FASTA: {e}")

        QTimer.singleShot(0, do_export)
            
            
    ############## CHECKBOXES TEXTS ############################
    def update_atg_checkbox_text(self):
        sender = self.sender()

        if sender == self.protein_checkbox_obj and self.protein_checkbox_obj.isChecked():
            self.protein_checkbox_obj.setEnabled(False)
            self.gene_checkbox_obj.setChecked(False)
            self.gene_checkbox_obj.setEnabled(True)
            self.atg_checkbox.setText("Check M")
        elif sender == self.gene_checkbox_obj and self.gene_checkbox_obj.isChecked():
            self.gene_checkbox_obj.setEnabled(False)
            self.protein_checkbox_obj.setChecked(False)
            self.protein_checkbox_obj.setEnabled(True)
            self.atg_checkbox.setText("Check ATG")
        elif not self.protein_checkbox_obj.isChecked() and not self.gene_checkbox_obj.isChecked():
            self.atg_checkbox.setText("Check")
            

    ############## RESET FILTERS AND SHOW ALL ############################
    def reset_filters(self):
        self.filter_textbox.clear()
        self.sequence_textbox.clear()
        self.accuracy_box.setValue(100)
        self.len_box_low.setValue(0)
        self.len_box_hi.setValue(900000)
        self.length_checkbox_obj.setChecked(False)
        self.protein_checkbox_obj.setChecked(True)
        self.gene_checkbox_obj.setChecked(False)
        self.atg_checkbox.setChecked(False)
        self.score_checkbox_obj.setChecked(False)
        self.eval_checkbox_obj.setChecked(False)

        self.genes_list.clear()
        self.len_box_low.setValue(self.min_seq_len)
        self.len_box_hi.setValue(self.max_seq_len)
        sequences = fetch_all_sequences("cleaned.db")
        for header, sequence in sequences:
            item = QListWidgetItem(f"{header} [{len(sequence)}]")
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            self.genes_list.addItem(item)

        self.amount_label.setText(f"Records: {self.genes_list.count()}")
        
        
    ############## HISTOGRAM CLICK ############################
    def connect_histogram_click(self):
        self.chart3.canvas.mpl_connect("pick_event", self.on_length_bar_click)
        
    ############## HISTOGRAM CLICK ############################
    def on_length_bar_click(self, event):
        if not hasattr(self, "length_bins") or not isinstance(event.artist, matplotlib.patches.Rectangle):
            return
        index = self.length_bar_patches.index(event.artist)
        bin_min = self.length_bins[index]
        bin_max = self.length_bins[index + 1] - 1
        self.len_box_low.setValue(bin_min)
        self.len_box_hi.setValue(bin_max)
    
    def view_selected_alignment(self):
        selected_items = self.query_list_widget.selectedItems()
        if not selected_items:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Warning", "Please select a query.")
            return
        query_name = selected_items[0].text()
        query_alignments = parse_blast_xml("output.xml")
        if query_name not in query_alignments:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Not Found", f"No alignments found for: {query_name}")
            return


        self.viewer_window = BlastAlignmentViewer(query_name, query_alignments[query_name])
        self.viewer_window.show()
        self.viewer_window.raise_()
        self.viewer_window.activateWindow()
        
    def info1(self):
        from PyQt6.QtWidgets import QMessageBox
        QMessageBox.information(self, "INFO", "To extract multiple queries, paste name fragments separated by ENTER (one per line). To search for multiple fragments within an individual name, paste the queries separated by commas.")
        
    ##################################################### BLAST UI ######################################################
    def setup_blast_tab(self):
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        left_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        right_layout = QVBoxLayout()
        right_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        outer_layout = QVBoxLayout()
        outer_layout.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        
        self.blast_title = QLabel("BLAST+")
        font_big = self.blast_title.font()
        font_big.setBold(True)
        font_big.setPointSize(20)
        self.blast_title.setFont(font_big)
        outer_layout.addWidget(self.blast_title, alignment=Qt.AlignmentFlag.AlignHCenter)
        outer_layout.addWidget(self.make_hor_separator())
        
        font_bold2 = self.label_analyse.font()
        font_bold2.setBold(True)
        font_bold2.setPointSize(12)
        #---------------------------LEFT---------------------------
        self.label_load = QLabel("Create BLAST+ database from FASTA file")
        self.label_load.setFont(font_bold2)
        
        load_base_layout = QHBoxLayout()
        load_base_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        self.load_base_button = QPushButton("LOAD CREATED DATABASE")
        self.load_base_button.setFixedSize(200, 40)
        self.load_base_button.clicked.connect(lambda: load_prev_database(self,self.label_database_name))
        
        self.open_blast_button = QPushButton("Import BLAST Result (XML)")
        self.open_blast_button.setFixedSize(200, 40)
        self.open_blast_button.clicked.connect(lambda: open_sequence_from_xml(self,self.rich_text_result))
        self.open_blast_button.setEnabled(False)
        
        database_layout = QHBoxLayout()
        database_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.base_select_button = QPushButton("Load Protein Database")
        self.base_select_button.setFixedSize(200, 40)
        self.base_select_button.clicked.connect(lambda: choose_database(self, self.label_database_name))
        self.base2_select_button = QPushButton("Load Gene Database")
        self.base2_select_button.setFixedSize(200, 40)
        self.base2_select_button.clicked.connect(lambda: choose_database_n(self, self.label_database_name))
        self.label_database_file = QLabel("Database Files:")
        self.label_database_file.setFont(font_bold2)
        self.label_database_name = QTextEdit("none")
        self.label_database_name.setReadOnly(True)
        self.label_database_name.setFixedSize(400,60)
        self.label_database_name.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.label_sequences = QLabel("Query sequences")
        self.label_sequences.setFont(font_bold2)
        self.rich_text_sequences = QTextEdit()
        self.rich_text_sequences.setPlaceholderText("Paste query sequences here")
        self.label_result = QLabel("Result")
        self.label_result.setFont(font_bold2)
        self.rich_text_result = QTextEdit()
        self.rich_text_sequences.setMinimumWidth(800)
        self.rich_text_result.setMinimumWidth(800)
        self.rich_text_sequences.setMinimumHeight(150)
        self.rich_text_result.setMinimumHeight(150)


        runblast_layout = QHBoxLayout()
        runblast_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.blast_button = QPushButton("Run BLASTP")
        self.blast_button.setFixedSize(150, 40)
        self.blast_button.clicked.connect(lambda: run_blast(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.blast_button.setEnabled(False)
        
        self.blastn_button = QPushButton("Run BLASTN")
        self.blastn_button.setFixedSize(150, 40)
        self.blastn_button.clicked.connect(lambda: run_blastn(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.blastn_button.setEnabled(False)
        
        self.blastx_button = QPushButton("Run BLASTX")
        self.blastx_button.setFixedSize(150, 40)
        self.blastx_button.clicked.connect(lambda: run_blastx(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.blastx_button.setEnabled(False)
        
        self.tblastn_button = QPushButton("Run TBLASTN")
        self.tblastn_button.setFixedSize(150, 40)
        self.tblastn_button.clicked.connect(lambda: run_tblastn(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.tblastn_button.setEnabled(False)
        
        self.tblastx_button = QPushButton("Run TBLASTX")
        self.tblastx_button.setFixedSize(150, 40)
        self.tblastx_button.clicked.connect(lambda: run_tblastx(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.tblastx_button.setEnabled(False)
        
        open_blast_layout = QHBoxLayout()
        open_blast_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        self.open_blast_file_button = QPushButton("Show Alignment")
        self.open_blast_file_button.setFixedSize(150, 40)
        self.open_blast_file_button.clicked.connect(lambda: show_blast_file("output.txt"))
        self.open_blast_file_button.setEnabled(False)
        
        self.save_sequences_button = QPushButton("Save Sequences")
        self.save_sequences_button.setFixedSize(150, 40)
        self.save_sequences_button.clicked.connect(lambda: save_sequences_from_blast(self))
        self.save_sequences_button.setEnabled(False)
        
        self.save_blast_xml_button = QPushButton("Save BLAST XML")
        self.save_blast_xml_button.setFixedSize(150, 40)
        self.save_blast_xml_button.clicked.connect(lambda:save_blast_output_as(self))
        self.save_blast_xml_button.setEnabled(False)
        
        #load_base_layout.addWidget(self.load_base_button)
        load_base_layout.addWidget(self.base_select_button)
        load_base_layout.addWidget(self.base2_select_button)
        open_blast_layout.addWidget(self.open_blast_button)

        runblast_layout.addWidget(self.blast_button)
        runblast_layout.addWidget(self.blastn_button)
        runblast_layout.addWidget(self.blastx_button)
        runblast_layout.addWidget(self.tblastn_button)
        runblast_layout.addWidget(self.tblastx_button)
        
        open_result_layout = QHBoxLayout()
        open_result_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        open_result_layout.addWidget(self.open_blast_file_button)
        open_result_layout.addWidget(self.save_sequences_button)
        open_result_layout.addWidget(self.save_blast_xml_button)
        
        widgets = [
            self.label_database_file,
            self.label_database_name, self.label_sequences,
            self.rich_text_sequences,
            
        ]
        
        left_layout.addWidget(self.label_load,alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addLayout(load_base_layout)
        left_layout.addLayout(database_layout)
        left_layout.addLayout(open_blast_layout)
        for widget in widgets:
            left_layout.addWidget(widget, alignment=Qt.AlignmentFlag.AlignHCenter)
        
        left_layout.addLayout(runblast_layout)
        
        widgets = [
            
            self.label_result, self.rich_text_result,
        ]
        
        for widget in widgets:
            left_layout.addWidget(widget, alignment=Qt.AlignmentFlag.AlignHCenter)
        
        left_layout.addLayout(open_result_layout)
        
        
        self.query_list_widget = QListWidget()
        self.query_list_widget.setMinimumHeight(150)
        self.view_alignment_button = QPushButton("View")
        self.view_alignment_button.setFixedSize(150, 40)
        self.view_alignment_button.clicked.connect(self.view_selected_alignment)
        right_layout.addWidget(QLabel("Query list"))
        right_layout.addWidget(self.query_list_widget)
        blast_info_layout = QHBoxLayout()
        blast_info_layout.addWidget(self.view_alignment_button)
        blast_info_layout.addStretch()
        blast_ref_layout = QVBoxLayout()
        blast_ref_layout.setAlignment(Qt.AlignmentFlag.AlignRight)
        ref1 = QLabel("BLAST+ tools are developed by the National Center for Biotechnology Information (NCBI)")
        ref2 = QLabel("(Camacho, 2009; Sayers, 2025)")
        ref1.setAlignment(Qt.AlignmentFlag.AlignRight)
        ref2.setAlignment(Qt.AlignmentFlag.AlignRight)
        blast_ref_layout.addWidget(ref1)
        blast_ref_layout.addWidget(ref2)
        
        
        blast_info_layout.addLayout(blast_ref_layout)
        
        
        right_layout.addLayout(blast_info_layout)
        
        #---------------------------RIGHT---------------------------
        #main_layout.addWidget(self.make_separator())
        main_layout.addLayout(left_layout, 5)
        main_layout.addWidget(self.make_separator())
        main_layout.addLayout(right_layout, 5)
        #main_layout.addWidget(self.make_separator())
        
        
        # --------------------------- QUERY LIST + VIEW BUTTON ---------------------------
        outer_layout.addLayout(main_layout)
        self.blast_tab.setLayout(outer_layout)



    
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    #window.resize(1600, 900)
    #self.setFixedSize(1600, 900)
    window.show()
    #window.showMaximized()
    sys.exit(app.exec())


    

    
    