from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
    QLineEdit, QTextEdit, QPushButton, QButtonGroup, QSpacerItem, QSizePolicy, QGraphicsView, QGraphicsScene, QGraphicsPixmapItem)
from PySide6.QtGui import QIcon, QPixmap, QPainter, QFont, QColor
from PySide6.QtCore import QProcess, Qt

import sys
import pandas as pd
import os
from operate_analysis import analyze

class ExerciseGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CurlyArrows")
        # Increase the window size to accommodate larger images
        self.setGeometry(100, 100, 900, 550)  
        # Updated icon path to use the GUI_images folder
        self.setWindowIcon(QIcon("GUI_images/icon_GUI.png"))
        self.selected_exercise = None  # Track selected exercise
        self.num_exercises = self.load_exercise_count()
        
        # Set Consolas font for the entire application
        self.font = QFont("Consolas")
        QApplication.instance().setFont(self.font)
        
        self.initUI()
        self.setStyle()
        # Display instructions on startup
        self.display_instructions()

    def load_exercise_count(self):
        try:
            df = pd.read_csv("exercise_model_answers.csv", sep=';')
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

        # Image display - increased by 1.5x
        self.image_label = QLabel("No Image Selected")
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_view = QGraphicsView()
        self.image_scene = QGraphicsScene()
        self.image_view.setScene(self.image_scene)
        
        # Increase the size by 1.5x
        self.image_view.setFixedSize(900, 250)
        
        # Disable scrollbars
        self.image_view.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.image_view.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        # Set rendering quality to high (fixed to use only available options)
        self.image_view.setRenderHint(QPainter.Antialiasing, True)
        self.image_view.setRenderHint(QPainter.SmoothPixmapTransform, True)
        
        # Set the background to match the application background
        self.image_view.setBackgroundBrush(QColor("#19232D"))
        
        top_layout.addWidget(self.image_view)

        # Input fields
        self.key_input = QLineEdit()
        self.key_input.setPlaceholderText("Enter OpenAI key")
        self.key_input.setFont(self.font)
        self.json_input = QLineEdit()
        self.json_input.setPlaceholderText("Enter JSON code")
        self.json_input.setFont(self.font)

        key_label = QLabel("OpenAI key")
        key_label.setFont(self.font)
        json_label = QLabel("JSON code")
        json_label.setFont(self.font)
        
        left_layout.addWidget(key_label)
        left_layout.addWidget(self.key_input)
        left_layout.addWidget(json_label)
        left_layout.addWidget(self.json_input)

        # Output text area
        self.output_text = QTextEdit()
        self.output_text.setPlaceholderText("Output will appear here...")
        self.output_text.setFont(self.font)
        left_layout.addWidget(self.output_text, 1)

        # Instructions button
        self.instructions_button = QPushButton("INSTRUCTIONS")
        self.instructions_button.setFixedSize(120, 30)
        self.instructions_button.setFont(self.font)
        self.instructions_button.clicked.connect(self.display_instructions)
        button_layout.addWidget(self.instructions_button)
        
        # Exercise label and buttons
        exercise_label = QLabel("Exercise")
        exercise_label.setFont(self.font)
        button_layout.addWidget(exercise_label)
        
        self.button_group = QButtonGroup()
        for i in range(1, self.num_exercises + 1):
            btn = QPushButton(str(i))
            btn.setFixedSize(40, 40)
            btn.setFont(self.font)
            self.button_group.addButton(btn, i)
            button_layout.addWidget(btn)
            btn.clicked.connect(self.set_selected_exercise)

        self.button_group.buttonClicked.connect(self.highlight_button)

        # Spacer
        button_layout.addItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Analyze button
        self.analyze_button = QPushButton("ANALYZE")
        self.analyze_button.setFixedSize(80, 30)
        self.analyze_button.setFont(self.font)
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
                background-color: #19232D;
                color: #F2F2F2;
                font-family: Consolas;
            }
            QLineEdit, QTextEdit {
                background-color: #222B35;
                border: 1px solid #F2F2F2;
                padding: 5px;
                font-family: Consolas;
                color: #A6A6A6;
            }
            QLineEdit[text=""], QTextEdit[plainText=""] {
                color: #A6A6A6;
            }
            QLineEdit:focus, QTextEdit:focus {
                color: #B0E686;
            }
            QLineEdit::placeholder, QTextEdit::placeholder {
                color: #A6A6A6;
            }
            QPushButton {
                background-color: #19232D;
                color: #F2F2F2;
                border: 1px solid #F2F2F2;
                font-family: Consolas;
            }
            QPushButton:checked {
                background-color: #ff4c4c;
                color: #F2F2F2;
            }
            QGraphicsView {
                border: none;
                background-color: #19232D;
            }
            #analyze_button {
                background-color: #4caf50;
                color: white;
                font-weight: bold;
            }
        """)
        
        # Set the analyze button style separately 
        self.analyze_button.setObjectName("analyze_button")

        # Update text colors programmatically to ensure they're applied
        self.key_input.setStyleSheet("color: #A6A6A6;")
        self.json_input.setStyleSheet("color: #A6A6A6;")
        self.output_text.setStyleSheet("color: #A6A6A6;")
        
        # Connect text changed signals to update color
        self.key_input.textChanged.connect(self.updateTextColor)
        self.json_input.textChanged.connect(self.updateTextColor)

    def updateTextColor(self):
        # Change text color to green when text is entered
        sender = self.sender()
        if sender.text():
            sender.setStyleSheet("color: #B0E686;")
        else:
            sender.setStyleSheet("color: #A6A6A6;")

    def set_selected_exercise(self):
        clicked_button = self.sender()
        self.selected_exercise = int(clicked_button.text())
        self.load_exercise_image()

    def display_instructions(self):
        # Clear any selected exercise button highlight
        for btn in self.button_group.buttons():
            btn.setStyleSheet("background-color: #19232D; color: #F2F2F2;")
        
        # Display instructions image
        image_path = "GUI_images/instructions.png"
        if os.path.exists(image_path):
            pixmap = QPixmap(image_path)
            self.image_scene.clear()
            # Add the image without scaling initially
            pixmap_item = self.image_scene.addPixmap(pixmap)
            # Reset transforms
            self.image_view.resetTransform()
            # Fit to view while preserving aspect ratio
            self.image_view.fitInView(pixmap_item, Qt.KeepAspectRatio)
        else:
            self.output_text.setText("Instructions image not found")

    def load_exercise_image(self):
        if self.selected_exercise is not None:
            image_path = f"GUI_images/exercise_{self.selected_exercise}.png"
            if os.path.exists(image_path):
                # Load image with highest quality
                pixmap = QPixmap(image_path)
                self.image_scene.clear()
                # Add the image to the scene
                pixmap_item = self.image_scene.addPixmap(pixmap)
                # Reset transforms
                self.image_view.resetTransform()
                # Fit to view while preserving aspect ratio
                self.image_view.fitInView(pixmap_item, Qt.KeepAspectRatio)
            else:
                self.output_text.setText("Image not found")

    def highlight_button(self, button):
        for btn in self.button_group.buttons():
            btn.setStyleSheet("background-color: #19232D; color: #F2F2F2; font-family: Consolas;")
        button.setStyleSheet("background-color: #ff4c4c; color: #F2F2F2; font-family: Consolas;")

    def run_analysis(self):
        ai_key = self.key_input.text() or False
        json_code = self.json_input.text()
        if not json_code or self.selected_exercise is None:
            self.output_text.setText("Error: Please enter JSON code and select an exercise.")
            return
        try:
            output = analyze(ai_key, self.selected_exercise, json_code)
            formatted_output = f"<p style='font-size:10pt; color: #B0E686; font-family: Consolas;'>{output}</p>"
            self.output_text.setHtml(formatted_output)
        except Exception as e:
            self.output_text.setText(f"Error: {str(e)}")

    def resizeEvent(self, event):
        # When window is resized, update the image view to maintain quality
        if hasattr(self, 'image_scene') and self.image_scene.items():
            pixmap_item = next((item for item in self.image_scene.items() if isinstance(item, QGraphicsPixmapItem)), None)
            if pixmap_item:
                self.image_view.fitInView(pixmap_item, Qt.KeepAspectRatio)
        super().resizeEvent(event)

    def closeEvent(self, event):
        self.deleteLater()
        event.accept()

if __name__ == "__main__":
    if QApplication.instance():
        app = QApplication.instance()
    else:
        app = QApplication(sys.argv)
    
    # Apply Consolas font to the entire application
    app.setFont(QFont("Consolas"))
    
    # For Windows taskbar icon to work properly
    # This needs to be done before creating the window
    if sys.platform == 'win32':
        import ctypes
        myappid = 'mycompany.curlyarrows.version0.1'  # Arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app.setWindowIcon(QIcon("GUI_images/icon_GUI.png"))
    window = ExerciseGUI()
    window.show()
    sys.exit(app.exec())