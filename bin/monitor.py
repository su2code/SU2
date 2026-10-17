# monitor: displays and updates the residuals of the su2 history file in a window.
import csv
import io
import os.path
import re
import sys

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5 import QtCore
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication, QCheckBox, QFileDialog, QFrame, QGridLayout, QGroupBox, \
                            QHBoxLayout, QLabel, QLineEdit, QMainWindow, QPushButton, QSplitter, \
                            QTabWidget, QVBoxLayout, QWidget

matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelsize'] = 20
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['axes.linewidth'] = 1.6
for axis in ('xtick', 'ytick'):
    matplotlib.rcParams[axis + '.major.size'] = 10.5
    matplotlib.rcParams[axis + '.major.width'] = 1.6
    matplotlib.rcParams[axis + '.minor.size'] = 6
    matplotlib.rcParams[axis + '.minor.width'] = 1.2
matplotlib.rcParams['legend.fontsize'] = 20

RESIDUAL_PATTERN = re.compile('rms|Res')
COEFFICIENT_PATTERN = re.compile('^(CD|CL|CMx|CMy|CMz|CFx|CFy|CFz)$')
MAX_POINTS = 4000
READ_BLOCK_BYTES = 32 * 1024 * 1024
UPDATE_INTERVAL_MS = 1000
LOGO_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logoSU2.png')


def cleanName(name):
    return name.replace('VARIABLES=', '').replace('"', '').replace(' ', '')


# ############################################################### #
# ############################################################### #
# ############################################################### #
class HistoryReader:
    """ incremental reader of the columns of an su2 history file matching a name pattern
    """

    def __init__(self, filename, pattern=RESIDUAL_PATTERN):
        self.filename = filename
        self.pattern = pattern
        self.reset()

    def reset(self):
        self.offset = 0
        self.cols = None
        self.names = []
        self.buf = np.empty((1024, 0))
        self.n = 0
        self.version = 0
        self.readHeader()

    # legacy files start with a TITLE line and have a ZONE line after the header
    def readHeader(self):
        try:
            with open(self.filename, 'rb') as f:
                first = f.readline()
                legacy = first.startswith(b'TITLE')
                header = f.readline() if legacy else first
                if legacy:
                    f.readline()
                if not header.endswith(b'\n'):
                    return
                offset = f.tell()
            allnames = next(csv.reader([header.decode()]))
        except (OSError, UnicodeDecodeError, csv.Error, StopIteration):
            return

        self.cols = [i for i, c in enumerate(allnames) if self.pattern.search(cleanName(c))]
        self.names = [cleanName(allnames[i]) for i in self.cols]
        self.buf = np.empty((1024, len(self.cols)))
        self.offset = offset

    # reads the complete lines appended since the last call; returns True if rows were added
    def update(self):
        try:
            size = os.path.getsize(self.filename)
        except OSError:
            return False

        if self.cols is not None and size < self.offset:
            self.reset()
        if self.cols is None:
            self.readHeader()
            if self.cols is None:
                return False

        added = False
        while size > self.offset:
            with open(self.filename, 'rb') as f:
                f.seek(self.offset)
                chunk = f.read(min(size - self.offset, READ_BLOCK_BYTES))
            end = chunk.rfind(b'\n')
            if end < 0:
                break
            try:
                rows = pd.read_csv(io.BytesIO(chunk[:end + 1]), header=None, usecols=self.cols,
                                   dtype=np.float64, on_bad_lines='skip').to_numpy()
            except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError):
                break
            self.offset += end + 1
            self.append(rows)
            added = True

        if added:
            self.version += 1
        return added

    def append(self, rows):
        needed = self.n + len(rows)
        if needed > len(self.buf):
            grown = np.empty((max(needed, 2 * len(self.buf)), self.buf.shape[1]))
            grown[:self.n] = self.buf[:self.n]
            self.buf = grown
        self.buf[self.n:needed] = rows
        self.n = needed


# ############################################################### #
# ############################################################### #
# ############################################################### #
class MatplotlibFigure(QWidget):
    """ setup of layout with canvas and toolbar for matplotlib figure
    """

    def __init__(self, pattern, ylabel, title):
        super().__init__()

        self.pattern = pattern
        self.ylabel = ylabel
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.reader = None
        self.lines = []
        self.visible = None
        self.lastkey = None

        layout = QGridLayout()
        self.setLayout(layout)

        groupbox = QGroupBox(title)
        layout.addWidget(groupbox)

        mpllayout = QVBoxLayout()
        groupbox.setLayout(mpllayout)
        mpllayout.addWidget(self.toolbar)
        mpllayout.addWidget(self.canvas)

    def setFile(self, filename):
        self.reader = HistoryReader(filename, self.pattern) if filename != "" else None
        self.visible = None
        self.lastkey = None
        self.figure.clf()
        self.lines = []
        self.canvas.draw_idle()

    def names(self):
        return self.reader.names if self.reader is not None else []

    # plots the residuals from iteration 'start' on; ymin/ymax of None means automatic range
    def plot(self, listonoff, normalize, start, ymin, ymax, update=True):
        if self.reader is None:
            return
        if update:
            self.reader.update()

        key = (tuple(listonoff), normalize, start, ymin, ymax, self.reader.version)
        if key == self.lastkey:
            return
        self.lastkey = key

        names = self.reader.names
        shown = [l for l in range(min(len(listonoff), len(names))) if listonoff[l]]
        ax = self.figure.axes[0] if self.figure.axes else self.figure.add_subplot(111)

        rebuild = self.visible != [names[l] for l in shown]
        if rebuild:
            self.visible = [names[l] for l in shown]
            ax.clear()
            ax.set_xlabel('iterations')
            ax.set_ylabel(self.ylabel)
            self.lines = [ax.plot([], [], label=names[l])[0] for l in shown]

        n = self.reader.n
        first = min(max(start, 1), max(n, 1))
        step = max(1, -(-(n - first) // MAX_POINTS))
        x = np.arange(first, n, step)
        for line, l in zip(self.lines, shown):
            y = self.reader.buf[first:n:step, l]
            if normalize and n > 1 and abs(self.reader.buf[1, l]) > 1e-6:
                y = 1 + y - self.reader.buf[1, l]
            line.set_data(x, y)

        ax.set_autoscaley_on(True)
        ax.relim()
        ax.autoscale_view()
        if ymin is not None or ymax is not None:
            low, high = ax.get_ylim()
            low = low if ymin is None else ymin
            high = high if ymax is None else ymax
            if low < high:
                ax.set_ylim(low, high)

        if rebuild:
            if self.lines:
                ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', prop={'size': 6, 'weight': 'bold'},
                          borderaxespad=None)
            self.figure.tight_layout()

        self.canvas.draw_idle()


# ############################################################### #
# ############################################################### #
# ############################################################### #
class MainApp(QMainWindow):
    """ Create the main window
    """
    def __init__(self):

        super().__init__()
        self.title = 'SU2 convergence monitor'

        # start at top left corner of the screen
        self.left = 0
        self.top = 0

        # set initial window height and width
        self.width = 1280
        self.height = 780

        # set the styles of fonts etc.
        self.setStyleSheet('font-size: 12pt')

        # user interface
        self.makeUI()

    def makeUI(self):
        # set window geometry
        self.setWindowTitle(self.title)
        self.setWindowIcon(QIcon(LOGO_FILE))
        self.setGeometry(self.left, self.top, self.width, self.height)

        # generate the tab layout
        self.table_widget = MyTableWidget()
        self.setCentralWidget(self.table_widget)

        self.show()


# ############################################################### #
# ############################################################### #
# ############################################################### #
class PlotTab(QWidget):
    """ one tab with the buttons, the data selection and the plot of a group of history columns
    """

    def __init__(self, owner, pattern, ylabel, title):
        super().__init__()
        self.owner = owner
        self.pattern = pattern
        self.ylabel = ylabel
        self.normalize = False
        self.chkbx = []
        self.chkbx_names = []
        self.make_tab(title)

    def setFile(self, filename, label, tooltip):
        self.show_filename.setText(label)
        self.show_filename.setToolTip(tooltip)
        self.figure.setFile(filename)
        self.update_chkbx(self.figure.names())
        self.plot_data()

    def setRunning(self, running):
        self.btn_plot_data.setText('stop' if running else 'start')
        self.btn_plot_data.setToolTip('Click to stop real time data gathering' if running
                                      else 'Click to start real time data gathering')

    # ################################### #
    # checkboxes for all the individual residual lines of the variables
    def update_chkbx(self, names):
        while self.chk_layout.count():
            widget = self.chk_layout.takeAt(0).widget()
            widget.setParent(None)
            widget.deleteLater()

        self.chkbx = []
        for name in names:
            box = QCheckBox(name)
            box.setChecked(True)
            box.toggled.connect(self.changedata)
            self.chk_layout.addWidget(box)
            self.chkbx.append(box)
        self.chkbx_names = list(names)

    # ################################### #
    def update_normalize(self, state):
        self.normalize = (state == QtCore.Qt.Checked)
        self.changedata()

    # ################################### #
    def make_tab(self, title):

        self.figure = MatplotlibFigure(self.pattern, self.ylabel, title)

        '''file dialog section'''
        self.btn_openfile = QPushButton("Select File")
        self.btn_openfile.clicked.connect(self.owner.getfile)
        self.show_filename = QLabel("No file selected")
        self.btn_openfile.setToolTip('click to select the file containing the residual data')

        '''button section'''
        # temporarily stop the realtime update
        self.btn_plot_data = QPushButton('stop')
        self.btn_plot_data.setToolTip('Click to stop real time data gathering')
        self.btn_plot_data.clicked.connect(self.owner.timer_startstop)

        self.btn_normalize = QCheckBox('Normalize')
        self.btn_normalize.setToolTip('Normalize all residuals with the value at iteration 1 ')
        self.btn_normalize.stateChanged.connect(self.update_normalize)

        '''plot range section'''
        self.edit_start = QLineEdit()
        self.edit_start.setPlaceholderText('1')
        self.edit_start.setToolTip('First iteration shown in the plot')
        self.edit_ymin = QLineEdit()
        self.edit_ymin.setPlaceholderText('auto')
        self.edit_ymin.setToolTip('Lower limit of the y-axis, empty for automatic range')
        self.edit_ymax = QLineEdit()
        self.edit_ymax.setPlaceholderText('auto')
        self.edit_ymax.setToolTip('Upper limit of the y-axis, empty for automatic range')

        range_layout = QGridLayout()
        for row, (label, edit) in enumerate([('start it.', self.edit_start),
                                             ('y min', self.edit_ymin),
                                             ('y max', self.edit_ymax)]):
            edit.editingFinished.connect(self.changedata)
            range_layout.addWidget(QLabel(label), row, 0)
            range_layout.addWidget(edit, row, 1)

        self.chk_layout = QVBoxLayout()

        # now we add all the buttons in a vertical list
        self.btn_layout = QVBoxLayout()
        self.btn_layout.addWidget(self.btn_openfile)
        self.btn_layout.addWidget(self.show_filename)
        self.btn_layout.addWidget(self.btn_plot_data)
        self.btn_layout.addWidget(self.btn_normalize)
        self.btn_layout.addLayout(range_layout)
        self.btn_layout.addLayout(self.chk_layout)

        # this adds stretch, so the buttons are always top-aligned (in a vbox)
        self.btn_layout.addStretch()

        '''left-side layout'''
        left = QFrame()
        left.setLayout(self.btn_layout)

        # combine the buttons and the canvas in a splitter layout
        splitter1 = QSplitter(QtCore.Qt.Horizontal)
        splitter1.addWidget(left)
        splitter1.addWidget(self.figure)
        splitter1.setSizes([75, 700])

        '''master layout of tab'''
        self.tab = QWidget()
        self.tab.layout = QHBoxLayout()
        self.tab.layout.addWidget(splitter1)
        self.tab.setLayout(self.tab.layout)

        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.tab)

    # ################################### #
    @staticmethod
    def parse(edit, convert):
        try:
            return convert(edit.text().strip())
        except ValueError:
            return None

    def draw(self, update):
        start = self.parse(self.edit_start, lambda t: int(float(t)))
        self.figure.plot([c.isChecked() for c in self.chkbx], self.normalize,
                         1 if start is None else start,
                         self.parse(self.edit_ymin, float), self.parse(self.edit_ymax, float),
                         update=update)

    # ################################### #
    def plot_data(self):
        if self.figure.reader is not None:
            self.figure.reader.update()
            if self.figure.names() != self.chkbx_names:
                self.update_chkbx(self.figure.names())
        self.draw(update=False)

    # ################################### #
    # redraw from the cached data after a change of the plot settings
    def changedata(self):
        self.draw(update=False)


# ############################################################### #
# ############################################################### #
# ############################################################### #
class MyTableWidget(QWidget):
    """ holds the tabs and the shared file name and update timer
    """

    def __init__(self):

        super().__init__()
        self.layout = QVBoxLayout(self)

        # initialization of the filename containing the history data
        self.filename = ""

        self.tabs = QTabWidget()
        self.tabs.resize(600, 600)

        self.plottabs = [PlotTab(self, RESIDUAL_PATTERN, 'residuals', 'Convergence'),
                         PlotTab(self, COEFFICIENT_PATTERN, 'coefficients', 'Coefficients')]
        self.tabs.addTab(self.plottabs[0], 'Convergence plot')
        self.tabs.addTab(self.plottabs[1], 'Coefficients plot')

        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        # setup a timer for plot update every 'UPDATE_INTERVAL_MS' ms
        self.timer = QtCore.QTimer()
        self.timer.setInterval(UPDATE_INTERVAL_MS)
        self.timer.timeout.connect(self.plot_data)
        self.timer.start()

    # ################################### #
    # file dialog for monitor input
    def getfile(self):

        fname = QFileDialog.getOpenFileName(self, 'Open file', 'history.csv', filter='csv files (*.csv)\n'
                                                                                    'dat files (*.dat)\n'
                                                                                    'all files (*)')
        if fname[0] == "":
            return

        self.filename = fname[0]
        fn = os.path.basename(self.filename)
        fb = os.path.basename(os.path.dirname(self.filename))
        for tab in self.plottabs:
            tab.setFile(self.filename, os.path.join(fb, fn), self.filename)

    # ################################### #
    # start or stop the realtime update of all tabs
    def timer_startstop(self):
        if self.timer.isActive():
            self.timer.stop()
        else:
            self.timer.start()
        for tab in self.plottabs:
            tab.setRunning(self.timer.isActive())

    # ################################### #
    def plot_data(self):
        for tab in self.plottabs:
            tab.plot_data()


# ############################################################### #
# ############################################################### #
# ############################################################### #
if __name__ == '__main__':

    def run_app():
        app = QApplication(sys.argv)
        app.setDesktopFileName('su2-monitor')
        app.setWindowIcon(QIcon(LOGO_FILE))
        mainWin = MainApp()
        mainWin.show()
        app.exec_()

    run_app()
