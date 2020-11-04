#pragma once

#include <Libs/Groom/Groom.h>
#include <QObject>

namespace shapeworks {

//! Qt Wrapper for Groom
/*!
 * The QGroom class wraps the Groom class to provide a QObject with a progress signal
 *
 */
class QGroom : public QObject, public Groom {

Q_OBJECT;

public:

  QGroom(ProjectHandle project);

protected:
  // override update_progress to emit q_signal
  void update_progress();

Q_SIGNALS:
  void progress(int);

};
}
