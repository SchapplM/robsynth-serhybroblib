% Jacobian time derivative of explicit kinematic constraints of
% fivebar1IC
% with respect to passive joint coordinates
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
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = fivebar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:08
% EndTime: 2020-04-27 06:19:08
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->4), mult. (8->6), div. (0->0), fcn. (4->4), ass. (0->5)
t20 = pkin(4) * (qJD(3) + qJD(4));
t19 = pkin(5) * (qJD(1) + qJD(2));
t18 = qJ(1) + qJ(2);
t17 = qJ(3) + qJ(4);
t1 = [cos(t18) * t19, cos(t17) * t20, 0; sin(t18) * t19, sin(t17) * t20, 0; 0, 0, 0;];
PhiD_p = t1;
