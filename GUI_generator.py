from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
    QLineEdit, QTextEdit, QPushButton, QButtonGroup, QSpacerItem, QSizePolicy, QGraphicsView, QGraphicsScene, QGraphicsPixmapItem)
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtCore import QProcess, Qt

import sys
import pandas as pd
import os
from operate_analysis import analyze

class ExerciseGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Curly Arrows 0.1")
        self.setGeometry(100, 100, 700, 400)
        self.setWindowIcon(QIcon("icon_GUI.png"))
        self.selected_exercise = None  # Track selected exercise
        self.num_exercises = self.load_exercise_count()
        self.initUI()
        self.setStyle()

    def load_exercise_count(self):
        try:
            df = pd.read_csv("exercise_model_answers.csv")
            return len(df)
        except Exception as e:
            print(f"Error loading CSV: {e}")
            return 9  # Default to 9 tiles if file fails to load

    def initUI(self):
        main_layout = QVBoxLayout(self)
        top_layout = QVBoxLayout()
        middle_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        button_layout = QVBoxLayout()

        # Image display
        self.image_label = QLabel("No Image Selected")
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_view = QGraphicsView()
        self.image_scene = QGraphicsScene()
        self.image_view.setScene(self.image_scene)
        self.image_view.setFixedSize(600, 300)  # Maintain aspect ratio
        top_layout.addWidget(self.image_view)

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

        # Exercise label and buttons
        button_layout.addWidget(QLabel("Exercise"))
        self.button_group = QButtonGroup()
        for i in range(1, self.num_exercises + 1):
            btn = QPushButton(str(i))
            btn.setFixedSize(40, 40)
            self.button_group.addButton(btn, i)
            button_layout.addWidget(btn)
            btn.clicked.connect(self.set_selected_exercise)

        self.button_group.buttonClicked.connect(self.highlight_button)

        # Spacer
        button_layout.addItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Analyze button
        self.analyze_button = QPushButton("ANALYZE")
        self.analyze_button.setFixedSize(80, 30)
        self.analyze_button.setStyleSheet("background-color: #4caf50; color: white; font-weight: bold;")
        self.analyze_button.clicked.connect(self.run_analysis)
        button_layout.addWidget(self.analyze_button)

        right_layout.addLayout(button_layout)
        middle_layout.addLayout(left_layout, 2)
        middle_layout.addLayout(right_layout, 1)

        main_layout.addLayout(top_layout)
        main_layout.addLayout(middle_layout)
        self.setLayout(main_layout)

    def setStyle(self):
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
        clicked_button = self.sender()
        self.selected_exercise = int(clicked_button.text())
        self.load_exercise_image()


    def load_exercise_image(self):
        if self.selected_exercise is not None:
            image_path = f"exercise_%s.png" %self.selected_exercise
            if os.path.exists(image_path):
                pixmap = QPixmap(image_path)
                self.image_scene.clear()
                self.image_scene.addPixmap(pixmap.scaled(600, 300, Qt.KeepAspectRatio,Qt.SmoothTransformation))
            else:
                self.output_text.setText("Image not found")

    def highlight_button(self, button):
        for btn in self.button_group.buttons():
            btn.setStyleSheet("background-color: #303030; color: white;")
        button.setStyleSheet("background-color: #ff4c4c; color: black;")

    def run_analysis(self):
        ai_key = self.key_input.text() or False
        json_code = self.json_input.text()
        if not json_code or self.selected_exercise is None:
            self.output_text.setText("Error: Please enter JSON code and select an exercise.")
            return
        try:
            output = analyze(ai_key, self.selected_exercise, json_code)
            formatted_output = f"<p style='font-size:10pt;'>{output}</p>"
            self.output_text.setHtml(formatted_output)
        except Exception as e:
            self.output_text.setText(f"Error: {str(e)}")

    def closeEvent(self, event):
        self.deleteLater()
        event.accept()

if __name__ == "__main__":
    if QApplication.instance():
        app = QApplication.instance()
    else:
        app = QApplication(sys.argv)
    window = ExerciseGUI()
    window.show()
    sys.exit(app.exec())
