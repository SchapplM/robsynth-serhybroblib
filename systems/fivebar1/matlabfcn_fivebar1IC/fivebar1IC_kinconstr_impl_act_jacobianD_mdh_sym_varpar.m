% Jacobian time derivative of explicit kinematic constraints of
% fivebar1IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% PhiD_a [(no of constraints)x(no. of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_a = fivebar1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:08
% EndTime: 2020-04-27 06:19:08
% DurationCPUTime: 0.05s
% Computational Cost: add. (12->8), mult. (16->12), div. (0->0), fcn. (8->8), ass. (0->7)
t12 = pkin(4) * (qJD(3) + qJD(4));
t11 = pkin(5) * (qJD(1) + qJD(2));
t10 = pkin(2) * qJD(1);
t9 = pkin(3) * qJD(3);
t8 = qJ(1) + qJ(2);
t7 = qJ(3) + qJ(4);
t1 = [cos(t8) * t11 - cos(qJ(1)) * t10, cos(t7) * t12 + cos(qJ(3)) * t9; sin(t8) * t11 - sin(qJ(1)) * t10, sin(t7) * t12 + sin(qJ(3)) * t9; 0, 0;];
PhiD_a = t1;
