% Jacobian time derivative of explicit kinematic constraints of
% fourbar1IC
% with respect to passive joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = fourbar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:48
% EndTime: 2020-04-24 20:15:48
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->6), mult. (8->6), div. (0->0), fcn. (4->4), ass. (0->4)
t11 = pkin(3) * (qJD(1) + qJD(2));
t10 = pkin(4) * qJD(3);
t9 = qJ(1) + qJ(2);
t1 = [-cos(t9) * t11, -cos(qJ(3)) * t10, 0; -sin(t9) * t11, -sin(qJ(3)) * t10, 0; 0, 0, 0;];
PhiD_p = t1;
