#include "GUIEdge.h"
#include "GUINode.h"
#include "GraphWidget.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

GUINode::GUINode(GraphWidget* graphWidget)
    : graph(graphWidget)
{
    //setFlag(ItemIsMovable); // A move of a GUINode should move the corresponding node.
    setFlag(ItemSendsGeometryChanges);
    setCacheMode(DeviceCoordinateCache);
    setZValue(-1);
}

void GUINode::addGUIEdge(GUIEdge* edge)
{
    edgeList << edge;
    edge->adjust();
}

QList<GUIEdge*> GUINode::edges() const
{
    return edgeList;
}
void GUINode::setID(size_t idx)
{
    this->idx=idx;
}

bool GUINode::advance()
{
    if (newPos == pos())
        return false;

    setPos(newPos);
    return true;
}

QRectF GUINode::boundingRect() const
{
    qreal adjust = 2;
    return QRectF( -15 - adjust, -15 - adjust, 33 + adjust, 33 + adjust);
}

QPainterPath GUINode::shape() const
{
    QPainterPath path;
    path.addEllipse(-15, -15, 30, 30);
    return path;
}

void GUINode::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget*)
{
    if(idx==1)
    {
        painter->setPen(Qt::NoPen);
      //    painter->setBrush(Qt::darkGray);
      //    painter->drawEllipse(-7, -7, 20, 20);

          QRadialGradient gradient(-3, -3, 10);
          if (option->state & QStyle::State_Sunken) {
              gradient.setCenter(3, 3);
              gradient.setFocalPoint(3, 3);
              gradient.setColorAt(1, QColor(Qt::red).light(120));
              gradient.setColorAt(0, QColor(Qt::red).light(120));
          } else {
              gradient.setColorAt(0, Qt::red);
              gradient.setColorAt(1, Qt::red);
          }
          painter->setBrush(gradient);

         // painter->setPen(QPen(Qt::black, 0));
          painter->drawEllipse(-15, -15, 30, 30);

    }
    else
    {
        painter->setPen(Qt::NoPen);
      //    painter->setBrush(Qt::darkGray);
      //    painter->drawEllipse(-7, -7, 20, 20);

          QRadialGradient gradient(-3, -3, 10);
          if (option->state & QStyle::State_Sunken) {
              gradient.setCenter(3, 3);
              gradient.setFocalPoint(3, 3);
              gradient.setColorAt(1, QColor(0,153,255).light(120));
              gradient.setColorAt(0, QColor(0,153,255).light(120));
          } else {
              gradient.setColorAt(0, QColor(0,153,255));
              gradient.setColorAt(1, QColor(0,153,255));
          }
          painter->setBrush(gradient);

         // painter->setPen(QPen(Qt::black, 0));
          painter->drawEllipse(-15, -15, 30, 30);
    }

    painter->setPen(Qt::white);
    painter->drawText(-6,8,QString::number(idx));
    setToolTip("Node: "+QString::number(idx));

}

QVariant GUINode::itemChange(GraphicsItemChange change, const QVariant& value)
{
    switch (change) {
    case ItemPositionHasChanged:
        foreach (GUIEdge* edge, edgeList)
            edge->adjust();
        break;
    default:
        break;
    };

    return QGraphicsItem::itemChange(change, value);
}
void GUINode::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
    ;
}
void GUINode::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    setToolTip(QString::number(idx));
    update();
    QGraphicsItem::mousePressEvent(event);
}

void GUINode::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    update();
    QGraphicsItem::mouseReleaseEvent(event);
}

