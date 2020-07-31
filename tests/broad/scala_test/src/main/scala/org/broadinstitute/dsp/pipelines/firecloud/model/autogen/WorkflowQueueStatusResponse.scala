/**
  * FireCloud
  * Genome analysis execution service.
  *
  * OpenAPI spec version: 0.1
  *
  *
  * NOTE: This class is auto generated by the swagger code generator program.
  * https://github.com/swagger-api/swagger-codegen.git
  * Do not edit the class manually.
  */
package org.broadinstitute.dsp.pipelines.firecloud.model.autogen

case class WorkflowQueueStatusResponse(
    /* estimated milliseconds until the current queue is submitted */
    estimatedQueueTimeMS: Integer,
    /* the number of workflows in the queue ahead of the user's first workflow */
    workflowsBeforeNextUserWorkflow: Integer,
    /* Map[String,Int] */
    workflowCountsByStatus: Any
)