#include "CompareWidget.h"

#include "ui_CompareWidget.h"

namespace shapeworks {

//---------------------------------------------------------------------------
CompareWidget::CompareWidget(QWidget *parent) : QWidget(parent), ui_(new Ui::CompareWidget) {
  ui_->setupUi(this);

  connect(ui_->surface_distance, &QRadioButton::toggled, this, &CompareWidget::settings_changed);
  connect(ui_->original, &QCheckBox::toggled, this, &CompareWidget::settings_changed);
  connect(ui_->groomed, &QCheckBox::toggled, this, &CompareWidget::settings_changed);
  connect(ui_->reconstructed, &QCheckBox::toggled, this, &CompareWidget::settings_changed);
  connect(ui_->original_slider, &QSlider::valueChanged, this, &CompareWidget::settings_changed);
  connect(ui_->groomed_slider, &QSlider::valueChanged, this, &CompareWidget::settings_changed);
  connect(ui_->reconstructed_slider, &QSlider::valueChanged, this, &CompareWidget::settings_changed);
}
//---------------------------------------------------------------------------

CompareWidget::~CompareWidget() { delete ui_; }

//---------------------------------------------------------------------------
CompareSettings CompareWidget::get_settings() {
  CompareSettings settings;
  settings.surface_distance_mode_ = ui_->surface_distance->isChecked();
  settings.original_checked_ = ui_->original->isChecked();
  settings.groomed_checked_ = ui_->groomed->isChecked();
  settings.reconstructed_checked_ = ui_->reconstructed->isChecked();
  settings.original_opacity_ = ui_->original_slider->value() / 100.0;
  settings.groomed_opacity_ = ui_->groomed_slider->value() / 100.0;
  settings.reconstructed_opacity_ = ui_->reconstructed_slider->value() / 100.0;
  return settings;
}

}  // namespace shapeworks
