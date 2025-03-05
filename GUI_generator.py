from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
    QLineEdit, QTextEdit, QPushButton, QButtonGroup, QSpacerItem, QSizePolicy)
from PySide6.QtGui import QIcon
from PySide6.QtCore import QProcess, Qt

import sys
from operate_analysis import analyze


class ExerciseGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Curly Arrows 0.1")
        self.setGeometry(100, 100, 700, 300)
        self.setWindowIcon(QIcon("icon_GUI.png"))
        self.selected_exercise = None  # Track selected exercise
        self.initUI()
        self.setStyle()

    def initUI(self):
        main_layout = QHBoxLayout(self)
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        button_layout = QVBoxLayout()  # Separate layout for buttons to align them properly

        # Instructions Button (Top Right)
        self.instructions_button = QPushButton("Instructions")
        self.instructions_button.clicked.connect(self.open_instructions)
        self.instructions_button.setFixedSize(100, 30)
        right_layout.addWidget(self.instructions_button, 0, Qt.AlignRight)  # Align to top right

        # Input fields
        self.key_input = QLineEdit()
        self.key_input.setPlaceholderText("Enter OpenAI key")
        self.json_input = QLineEdit()
        self.json_input.setPlaceholderText("Enter JSON code")

        left_layout.addWidget(QLabel("OpenAI key"))
        left_layout.addWidget(self.key_input)
        left_layout.addWidget(QLabel("JSON code"))
        left_layout.addWidget(self.json_input)

        # Output text area
        self.output_text = QTextEdit()
        self.output_text.setPlaceholderText("Output will appear here...")
        left_layout.addWidget(self.output_text, 1)

        # Exercise label and buttons in a vertical layout
        button_layout.addWidget(QLabel("Exercise"))
        self.button_group = QButtonGroup()
        for i in range(1, 10):  # Creates buttons 1-9
            btn = QPushButton(str(i))
            btn.setFixedSize(40, 40)  # Make buttons square
            self.button_group.addButton(btn, i)
            button_layout.addWidget(btn)
            btn.clicked.connect(self.set_selected_exercise)

        self.button_group.buttonClicked.connect(self.highlight_button)

        # Spacer to push analyze button to the bottom
        button_layout.addItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Analyze button
        self.analyze_button = QPushButton("ANALYZE")
        self.analyze_button.setFixedSize(80, 30)
        self.analyze_button.setStyleSheet("background-color: #4caf50; color: white; font-weight: bold;")
        self.analyze_button.clicked.connect(self.run_analysis)
        button_layout.addWidget(self.analyze_button)

        right_layout.addLayout(button_layout)  # Add button layout to right layout

        main_layout.addLayout(left_layout, 2)
        main_layout.addLayout(right_layout, 1)
        self.setLayout(main_layout)

    def setStyle(self):
        """Apply the dark theme and button styles"""
        self.setStyleSheet("""
            QWidget {
                background-color: #303030;
                color: white;
            }
            QLineEdit, QTextEdit {
                background-color: #404040;
                color: white;
                border: 1px solid gray;
                padding: 5px;
            }
            QPushButton {
                background-color: #303030;
                color: white;
                border: 1px solid gray;
            }
            QPushButton:checked {
                background-color: #ff4c4c;
                color: black;
            }
        """)

    def set_selected_exercise(self):
        """Store the selected exercise number"""
        clicked_button = self.sender()
        self.selected_exercise = int(clicked_button.text())

    def highlight_button(self, button):
        """Highlight the selected button"""
        for btn in self.button_group.buttons():
            btn.setStyleSheet("background-color: #303030; color: white;")
        button.setStyleSheet("background-color: #ff4c4c; color: black;")

    def run_analysis(self):
        """Run analysis after checking input validity"""
        ai_key = self.key_input.text() or False  # Set to False if empty
        json_code = self.json_input.text()
        
        if not json_code or self.selected_exercise is None:
            self.output_text.setText("Error: Please enter JSON code and select an exercise.")
            return

        output = analyze(ai_key, self.selected_exercise, json_code)
        formatted_output = f"<p style='font-size:10pt;'>{output}</p>"
        self.output_text.setHtml(formatted_output)

    def open_instructions(self):
        """Open the instructions PDF"""
        QProcess.startDetached("xdg-open", ["instruction_CA.pdf"])  # Linux
        # Uncomment for Windows/macOS
        # QProcess.startDetached("open", ["instruction_CA.pdf"])  # macOS
        # QProcess.startDetached("start", ["instruction_CA.pdf"])  # Windows

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ExerciseGUI()
    app.setWindowIcon(QIcon("icon_GUI.png"))
    window.show()
    sys.exit(app.exec())
